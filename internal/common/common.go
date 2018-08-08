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

package common

import (
	"encoding/binary"
	"fmt"
	"io"

	"github.com/googlegenomics/htsget/internal/genomics"
)

const (
	maximumReadLength = 1 << 29

	// This ID is used as a virtual bin ID for (unused) chunk metadata.
	MetadataBeanID = 37450
)

// ChecksMagic checks the magic bytes from the provided reader.
func CheckMagic(r io.Reader, want []byte) error {
	got := make([]byte, len(want))
	if err := Read(r, &got); err != nil {
		return fmt.Errorf("reading magic: %v", err)
	}
	for i, n := 0, len(want); i < n; i++ {
		if got[i] != want[i] {
			return fmt.Errorf("wrong magic %v (wanted %v)", got, want)
		}
	}
	return nil
}

// Read reads the value from the provided reader into the provided interface.
func Read(r io.Reader, v interface{}) error {
	return binary.Read(r, binary.LittleEndian, v)
}

// RegionContainsBin indicates if the given region contains the bin described by
// referenceID and binID.
func RegionContainsBin(region genomics.Region, referenceID int32, binID uint32, bins []uint16) bool {
	if region.ReferenceID >= 0 && referenceID != region.ReferenceID {
		return false
	}

	if region.Start == 0 && region.End == 0 {
		return true
	}

	for _, id := range bins {
		if uint32(id) == binID {
			return true
		}
	}
	return false
}

// BinsForRange calculates the list of bins that may overlap with region [beg,end) (zero-based).
// This function is derived from the C examples in the CSI index specification.
func BinsForRange(start, end uint32, minShift, depth int32) []uint16 {
	if end == 0 || end > maximumReadLength {
		end = maximumReadLength
	}
	if end <= start {
		return nil
	}
	if start > maximumReadLength {
		return nil
	}

	end--
	var bins []uint16
	for l, t, s := uint(0), uint(0), uint(minShift+depth*3); l <= uint(depth); l++ {
		b := t + (uint(start) >> s)
		e := t + (uint(end) >> s)
		for i := b; i <= e; i++ {
			bins = append(bins, uint16(i))
		}
		s -= 3
		t += 1 << (l * 3)
	}
	return bins
}
