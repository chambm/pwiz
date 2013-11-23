﻿/*
 * Original author: Nicholas Shulman <nicksh .at. u.washington.edu>,
 *                  MacCoss Lab, Department of Genome Sciences, UW
 *
 * Copyright 2012 University of Washington - Seattle, WA
 *
 * Licensed under the Apache License, Version 2.0 (the "License");
 * you may not use this file except in compliance with the License.
 * You may obtain a copy of the License at
 *
 *     http://www.apache.org/licenses/LICENSE-2.0
 *
 * Unless required by applicable law or agreed to in writing, software
 * distributed under the License is distributed on an "AS IS" BASIS,
 * WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
 * See the License for the specific language governing permissions and
 * limitations under the License.
 */

using System;
using System.Collections.Generic;
using System.ComponentModel;
using pwiz.Common.DataBinding.Attributes;
using pwiz.Skyline.Model.Databinding.Collections;
using pwiz.Skyline.Model.DocSettings;

namespace pwiz.Skyline.Model.Databinding.Entities
{
    // ReSharper disable LocalizableElement
    [AnnotationTarget(AnnotationDef.AnnotationTarget.protein)]
    public class Protein : SkylineDocNode<PeptideGroupDocNode>
    {
        public Protein(SkylineDataSchema dataSchema, IdentityPath identityPath) : base(dataSchema, identityPath)
        {
        }

        private Peptides _peptides;
        [OneToMany(ForeignKey = "Protein")]
        public Peptides Peptides
        {
            get { return _peptides = _peptides ?? new Peptides(this); }
        }

        [OneToMany(ItemDisplayName = "ResultFile")]
        [Obsolete]
        public IDictionary<ResultKey, ResultFile> Results
        {
            get
            {
                var dict = new Dictionary<ResultKey, ResultFile>();
                var document = SrmDocument;
                if (document.Settings.HasResults)
                {
                    for (int iResult = 0; iResult < document.Settings.MeasuredResults.Chromatograms.Count; iResult++)
                    {
                        var replicate = new Replicate(DataSchema, iResult);
                        var chromatogramSet = document.Settings.MeasuredResults.Chromatograms[iResult];
                        for (int iFile = 0; iFile < chromatogramSet.MSDataFileInfos.Count; iFile++)
                        {
                            var resultFile = new ResultFile(replicate, chromatogramSet.MSDataFileInfos[iFile].FileId, 0);
                            dict.Add(new ResultKey(replicate, iFile), resultFile);
                        }
                    }
                }
                return dict;
            }
        }

        protected override PeptideGroupDocNode CreateEmptyNode()
        {
            return new PeptideGroupDocNode(new PeptideGroup(), null, null, new PeptideDocNode[0]);
        }

        [DisplayName("ProteinName")]
        public string Name
        {
            get { return DocNode.Name; }
            set { ChangeDocNode(DocNode.ChangeName(value)); }
        }
        [DisplayName("ProteinDescription")]
        public string Description
        {
            get { return DocNode.Description; }
            set { ChangeDocNode(DocNode.ChangeDescription(value)); }
        }
        [DisplayName("ProteinSequence")]
        public string Sequence 
        {
            get { return DocNode.PeptideGroup.Sequence; }
        }
        [DisplayName("ProteinNote")]
        public string Note
        {
            get { return DocNode.Note; }
            set { ChangeDocNode(DocNode.ChangeAnnotations(DocNode.Annotations.ChangeNote(value)));}
        }
        public override string ToString()
        {
            return Name;
        }
    }
}
