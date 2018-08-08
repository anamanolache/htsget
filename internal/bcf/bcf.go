// Copyright 2018 Google Inc.
//
// Licensed under the Apache License, Version 2.0 (the "License");
// you may not use this file except in compliance with the License.
// You may obtain a copy of the License at
//
// https://www.apache.org/licenses/LICENSE-2.0
//
// Unless required by applicable law or agreed to in writing, software
// distributed under the License is distributed on an "AS IS" BASIS,
// WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
// See the License for the specific language governing permissions and
// limitations under the License.

package bcf

import (
	"bufio"
	"compress/gzip"
	"fmt"
	"io"
	"io/ioutil"
	"strconv"
	"strings"

	"github.com/googlegenomics/htsget/internal/bgzf"
	"github.com/googlegenomics/htsget/internal/binary"
	"github.com/googlegenomics/htsget/internal/csi"
	"github.com/googlegenomics/htsget/internal/genomics"
)

const (
	bcfMagic = "BCF\x02\x02"
	csiMagic = "CSI\x01"
)

// GetReferenceID retrieves the reference id of the given referenceName
// from the provided bcf file.
func GetReferenceID(bcf io.Reader, referenceName string) (int, error) {
	gzr, err := gzip.NewReader(bcf)
	if err != nil {
		return 0, fmt.Errorf("initializing gzip reader: %v", err)
	}
	defer gzr.Close()

	if err := binary.ExpectBytes(gzr, []byte(bcfMagic)); err != nil {
		return 0, fmt.Errorf("checking magic of BCF file: %v", err)
	}

	var length uint32
	if err := binary.Read(gzr, &length); err != nil {
		return 0, fmt.Errorf("reading header length: %v", err)
	}

	headerReader := io.LimitReader(gzr, int64(length))
	scanner := bufio.NewScanner(headerReader)
	var id int
	var contigsFound bool
	for scanner.Scan() {
		if line := scanner.Text(); strings.HasPrefix(line, "##contig") {
			contigsFound = true
			if contigDefinesReference(line, referenceName) {
				idx, err := getIdx(line)
				if err != nil {
					return 0, fmt.Errorf("getting idx: %v", err)
				}
				if idx > -1 {
					return idx, nil
				}
				return id, nil
			}
			id++
		} else {
			if contigsFound {
				return 0, fmt.Errorf("reference name not found")
			}
		}
	}
	if err := scanner.Err(); err != nil {
		return 0, fmt.Errorf("scanning header: %v", err)
	}
	return 0, fmt.Errorf("region id not found")
}

func contigDefinesReference(contig, refName string) bool {
	index := strings.Index(contig, fmt.Sprintf("ID=%s", refName))
	if index == -1 {
		return false
	}
	if nextChr := contig[index+len("ID=")+len(refName)]; nextChr != ',' && nextChr != '>' {
		return false
	}
	return true
}

func getIdx(contig string) (int, error) {
	index := strings.Index(contig, "IDX=")
	if index == -1 {
		return -1, nil
	}
	index += len("IDX=")
	var buff []byte
	for n := len(contig); index < n; index++ {
		chr := contig[index]
		if chr == ',' || chr == '>' {
			break
		}
		buff = append(buff, chr)
	}
	idx, err := strconv.Atoi(string(buff))
	if err != nil {
		return -1, fmt.Errorf("parsing IDX value: %v", err)
	}
	return idx, nil
}

// Read reads index data from csi and returns a set of BGZF chunks covering
// the header and all mapped reads that fall inside the specified region.  The
// first chunk is always the BCF header.
func Read(csiFile io.Reader, region genomics.Region) ([]*bgzf.Chunk, error) {
	gzr, err := gzip.NewReader(csiFile)
	if err != nil {
		return nil, fmt.Errorf("initializing gzip reader: %v", err)
	}
	defer gzr.Close()
	if err := binary.ExpectBytes(gzr, []byte(csiMagic)); err != nil {
		return nil, fmt.Errorf("checking magic: %v", err)
	}

	var minShift int32
	if err := binary.Read(gzr, &minShift); err != nil {
		return nil, fmt.Errorf("reading # bits for the minimal interval (min_shift): %v", err)
	}
	var depth int32
	if err := binary.Read(gzr, &depth); err != nil {
		return nil, fmt.Errorf("reading depth of binary index: %v", err)
	}
	bins := csi.BinsForRange(region.Start, region.End, minShift, depth)

	var laux int32
	if err := binary.Read(gzr, &laux); err != nil {
		return nil, fmt.Errorf("reading length of auxiliary data: %v", err)
	}
	if _, err := io.CopyN(ioutil.Discard, gzr, int64(laux)); err != nil {
		return nil, fmt.Errorf("reading past auxiliary data: %v", err)
	}

	header := &bgzf.Chunk{End: bgzf.LastAddress}
	chunks := []*bgzf.Chunk{header}
	var refCount int32
	if err := binary.Read(gzr, &refCount); err != nil {
		return nil, fmt.Errorf("reading the number of reference sequences: %v", err)
	}
	for reference := int32(0); reference < refCount; reference++ {
		var binCount int32
		if err := binary.Read(gzr, &binCount); err != nil {
			return nil, fmt.Errorf("reading bin count: %v", err)
		}
		for j := int32(0); j < binCount; j++ {
			var bin struct {
				ID     uint32
				Offset uint64
				Chunks int32
			}
			if err := binary.Read(gzr, &bin); err != nil {
				return nil, fmt.Errorf("reading bin header: %v", err)
			}

			includeChunks := csi.RegionContainsBin(region, reference, bin.ID, bins)
			for k := int32(0); k < bin.Chunks; k++ {
				var chunk bgzf.Chunk
				if err := binary.Read(gzr, &chunk); err != nil {
					return nil, fmt.Errorf("reading chunk: %v", err)
				}
				if bin.ID == csi.MetadataBeanID {
					continue
				}
				if includeChunks && (chunk.End >= bgzf.Address(bin.Offset)) {
					chunks = append(chunks, &chunk)
				}
				if header.End > chunk.Start {
					header.End = chunk.Start
				}
			}
		}
	}
	return chunks, nil
}
