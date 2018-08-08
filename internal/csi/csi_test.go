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

package csi

import (
	"fmt"
	"os"
	"testing"

	"github.com/googlegenomics/htsget/internal/genomics"
)

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
