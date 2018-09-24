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

package api

import (
	"context"
	"fmt"

	"io"

	"cloud.google.com/go/storage"
	"github.com/googlegenomics/htsget/internal/bgzf"
	"github.com/googlegenomics/htsget/internal/genomics"
)

type chunksRequest struct {
	indexObjects   []*storage.ObjectHandle
	blockSizeLimit uint64
	region         genomics.Region
	read           func(csiFile io.Reader, region genomics.Region) ([]*bgzf.Chunk, error)
}

func (req *chunksRequest) handle(ctx context.Context) ([]*bgzf.Chunk, error) {
	var index *storage.Reader
	var err error
	for _, object := range req.indexObjects {
		index, err = object.NewReader(ctx)
		if err == nil {
			break
		}
	}
	if err != nil {
		return nil, newStorageError("opening index", err)
	}
	defer index.Close()

	chunks, err := req.read(index, req.region)
	if err != nil {
		return nil, fmt.Errorf("reading index: %v", err)
	}
	return bgzf.Merge(chunks, req.blockSizeLimit), nil
}
