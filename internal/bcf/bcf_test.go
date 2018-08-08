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
	"fmt"
	"os"
	"strings"
	"testing"

	"github.com/googlegenomics/htsget/internal/genomics"
)

func TestGetReferenceId(t *testing.T) {
	testCases := []struct {
		file   string
		name   string
		id     int
		errMsg string
	}{
		// The test file bcf_with_idx.bcf.gz has the chr2 line before the chr1 line
		// but the region id is given by the IDX field.
		{"bcf_with_idx.bcf.gz", "chr1", 0, ""},
		{"bcf_with_idx.bcf.gz", "chr2", 1, ""},
		{"bcf_with_idx.bcf.gz", "chr10", 9, ""},
		{"bcf_without_idx.bcf.gz", "19", 0, ""},
		{"bcf_without_idx.bcf.gz", "Y", 2, ""},
		{"bcf_without_idx.bcf.gz", "Z", 0, "reference name not found"},
	}

	for _, tc := range testCases {
		t.Run(fmt.Sprintf("%s_%s", tc.file, tc.name), func(t *testing.T) {
			r, err := os.Open(fmt.Sprintf("testdata/%s", tc.file))
			if err != nil {
				t.Fatalf("Failed to open testdata: %v", err)
			}
			defer r.Close()

			if id, err := GetReferenceID(r, tc.name); err != nil && !strings.Contains(err.Error(), tc.errMsg) {
				t.Fatalf("GetReferenceID() returned unexpected error: %v", err)
			} else if id != tc.id {
				t.Fatalf("Wrong reference ID: got %d, want %d", id, tc.id)
			}
		})
	}
}

func TestContigDefinesReference(t *testing.T) {
	testCases := []struct {
		line string
		ref  string
		want bool
	}{
		{"##contig=<ID=chr1,length=248956422,IDX=0>", "chr1", true},
		{"##contig=<ID=chr10,length=248956422,IDX=0>", "chr1", false},
		{"##contig=<ID=Y,length=248956422,IDX=0>", "chr1", false},
		{"##contig=<length=248956422,IDX=0>", "Y", false},
	}

	for i, tc := range testCases {
		t.Run(string(i), func(t *testing.T) {
			if got := contigDefinesReference(tc.line, tc.ref); got != tc.want {
				t.Fatalf("Wrong contigDefinesReference response, want %v, got %v ", tc.want, got)
			}
		})
	}
}

func TestGetIDX(t *testing.T) {
	testCases := []struct {
		line string
		want int
	}{
		{"##contig=<ID=chr1,length=248956422>", -1},
		{"##contig=<ID=chr1,length=248956422,IDX=0>", 0},
		{"##contig=<ID=chr1,length=248956422,IDX=7>", 7},
		{"##contig=<ID=chr1,length=248956422,IDX=125>", 125},
		{"##contig=<ID=chr1,IDX=125,length=248956422>", 125},
	}

	for _, tc := range testCases {
		t.Run(tc.line, func(t *testing.T) {
			if got, _ := getIdx(tc.line); got != tc.want {
				t.Fatalf("Wrong getIdx response, want %d, got %d ", tc.want, got)
			}
		})
	}
}

func TestRegionRead(t *testing.T) {
	testCases := []struct {
		refId  int32
		start  uint32
		end    uint32
		chunks int
	}{
		{0, 93822816, 93825705, 2},
		{0, 0, 1000000, 1},
	}

	for _, tc := range testCases {
		t.Run(fmt.Sprintf("%d(%d-%d)", tc.refId, tc.start, tc.end), func(t *testing.T) {
			r, err := os.Open("testdata/sample.bcf.gz.csi")
			if err != nil {
				t.Fatalf("Failed to open testdata: %v", err)
			}
			defer r.Close()

			region := genomics.Region{
				ReferenceID: tc.refId,
				Start:       tc.start,
				End:         tc.end,
			}
			chunks, err := Read(r, region)
			if err != nil {
				t.Fatalf("Read() returned unexpected error: %v", err)
			}
			if got, want := len(chunks), tc.chunks; got != want {
				t.Fatalf("Wrong number of chunks: got %d, want %d", got, want)
			}
		})
	}
}
