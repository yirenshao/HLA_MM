workflow major_hla_mm {
    
    File host_hla_typing
    File donor_hla_typing
    File ref_hla_prot_fasta
    File hlathena_avail_alleles

    

    # Runtime parameters
    Int preemptible = 3 # non-negative interger value for preemptible 0 means not preemptible,
                        # otherwise 1,2,... is the max number of pre-emptible tries
    Int memoryGB = 8
    Int diskGB = 10

    # Workflow parameters
    String runID #sample_name
    #File alleles_file
    #Array[String] alleles = read_lines(alleles_file)

    #File peptide_list
    String peptide_col_name

    Boolean exists_ctex
    Boolean exists_expr
    String expr_col_name
    Boolean logtransform_expr
    Boolean aggregate_pep

    Array[String] lens
    String models_path = "gs://msmodels/"
    String features_all = "features_AAPos_AAPCA_LogTPM_CNN_Kidera_Gene"
    String feature_sets = "features_AAPos_AAPCA_Kidera"

    String assign_by_ranks_or_scores
    String assign_threshold
    String assign_colors

    call HLAthena_preprocessing {
        input:
            host_hla_typing = host_hla_typing,
            donor_hla_typing = donor_hla_typing,
            ref_hla_prot_fasta = ref_hla_prot_fasta,
            runID = runID,
            hlathena_avail_alleles = hlathena_avail_alleles

    }

    ### Parse HLA alleles
    call parse_alleles_task {
        input:
            preemptible = preemptible,
            memoryGB = memoryGB,
            diskGB = diskGB,
                    
            alleles_file = HLAthena_preprocessing.HLAthena_input_hla
    }
    Array[String] alleles_parsed = read_lines(parse_alleles_task.alleles_parsed)

    ### Check whether the input is in fasta format and convert to peptides accordingly
    call check_for_fasta_and_conver_task {
        input:
            preemptible = preemptible,
            memoryGB = memoryGB,
            diskGB = diskGB,

            peptide_list = HLAthena_preprocessing.HLAthena_input_peptides,
            peptide_col_name = peptide_col_name,
            lens = lens
    }

    ### If up/dn context is in the input file, compute and append cleavability score
    if (exists_ctex) {
        call predict_cleavability_task {
            input:
                preemptible = preemptible,
                memoryGB = memoryGB,
                diskGB = diskGB,

                peptide_list = check_for_fasta_and_conver_task.peptide_list_out,
                peptide_col_name = peptide_col_name
        }
    }

    ### Optionally aggregate peptides based on a given column (e.g. agregate all possible transcript precursors of a peptide)
    if (aggregate_pep) {
        call aggregate_pep_task {
            input:
                preemptible = preemptible,
                memoryGB = memoryGB,
                diskGB = diskGB,

                peptide_list = if exists_ctex then predict_cleavability_task.peptide_list_clev else check_for_fasta_and_conver_task.peptide_list_out,
                peptide_col_name = peptide_col_name,
                exists_expr = exists_expr,
                expr_col_name = expr_col_name
        }
    }

    ### If expression is in the input file, log-transform as needed
    File? peptide_list_cond = if (aggregate_pep) then aggregate_pep_task.peptide_list_agg else ( if exists_ctex then predict_cleavability_task.peptide_list_clev else check_for_fasta_and_conver_task.peptide_list_out)
    if (exists_expr) {
        call transform_expression_task {
            input:
                preemptible = preemptible,
                memoryGB = memoryGB,
                diskGB = diskGB,

                peptide_list = peptide_list_cond,
                expr_col_name = expr_col_name,
                logtransform_expr = logtransform_expr
        }
    }

    ### Split merged peptides by length
    call split_peptides_len_task {
        input:
            preemptible = preemptible,
            memoryGB = memoryGB,
            diskGB = diskGB,

            peptide_list = if (exists_expr) then transform_expression_task.peptide_list_clev_expr else peptide_list_cond,
            peptide_col_name = peptide_col_name,
            lens = lens
    }
    Array[String] lens_present = read_lines(split_peptides_len_task.lens_present)

    ### Generate peptide sequences features: dummy encoding only
    ### (blosum and fuzzy generated on demand from the dummy encoding upon prediction)
    scatter (pepfilelen in split_peptides_len_task.pep_len_files) {
        call featurize_encoding_task {
            input:
                preemptible = preemptible,
                memoryGB = memoryGB,
                diskGB = diskGB,

                peptide_len_list = pepfilelen,
                encoding = "dummy",
                peptide_col_name = peptide_col_name
        }
    }

    ### Predict with MS models (now includes ranks)
    scatter (allele_featfile in cross(alleles_parsed, featurize_encoding_task.feature_files)) {
        call predict_ms_allele_task {
            input:
                preemptible = preemptible,
                memoryGB = memoryGB,
                diskGB = diskGB,

                patient = runID,
                featfile = allele_featfile.right,
                lens = lens_present,
                exists_ctex = exists_ctex,
                exists_expr = exists_expr,
                allele = allele_featfile.left,
                models_path = models_path,
                peptide_col_name = peptide_col_name,
                features_all = features_all,
                feature_sets = feature_sets
        }
    }

    ### Merge predictions into ensemble model scores: MSIntrinsic, MSIntrinsicC, MSIntrinsicEC
    scatter (len in lens_present) {
        call merge_mspreds_task {
            input:
                preemptible = preemptible,
                memoryGB = memoryGB,
                diskGB = diskGB,

                patient = runID,
                alleles = alleles_parsed,
                len = len,
                exists_ctex = exists_ctex,
                exists_expr = exists_expr,
                mspred_files = predict_ms_allele_task.mspred_files,
                peptide_col_name = peptide_col_name
        }
    }

    ### Convert scores to ranks
    scatter (preds_file in merge_mspreds_task.mspreds_file_wide) {
        call get_ranks_task {
            input:
                preemptible = preemptible,
                memoryGB = memoryGB,
                diskGB = diskGB,

                patient = runID,
                alleles = alleles_parsed,
                lens = lens_present,
                exists_ctex = exists_ctex,
                exists_expr = exists_expr,
                models_path = models_path,
                peptide_col_name = peptide_col_name,
                preds_file = preds_file
        }
    }

    ### Concatenate ranks files for all lengths
    call ranks_concat_lens_task {
        input:
            preemptible = preemptible,
            memoryGB = memoryGB,
            diskGB = diskGB,

            patient = runID,
            peptide_col_name = peptide_col_name,
            exists_ctex = exists_ctex,
            exists_expr = exists_expr,
            assign_by_ranks_or_scores = assign_by_ranks_or_scores,
            assign_threshold = assign_threshold,
            assign_colors = assign_colors,
            ranks_len_files = get_ranks_task.ranks_file
    }


    ### Workflow level outputs
    output {
        ranks_concat_lens_task.sample_predictions
        ranks_concat_lens_task.sample_allele_assignment_counts
        ranks_concat_lens_task.sample_allele_assignment_counts_for_plot
        ranks_concat_lens_task.sample_plots
    }
}



task HLAthena_preprocessing{
    File host_hla_typing
    File donor_hla_typing
    File ref_hla_prot_fasta 
    File hlathena_avail_alleles
    String runID

    command {

                python3 /scripts/major_HLA_MM_HLAthena_preprocessing.py -donor ${donor_hla_typing} -host ${host_hla_typing} -ref ${ref_hla_prot_fasta} -HLAthena_avail_alleles ${hlathena_avail_alleles} -sample_name ${runID}
        }

        runtime {
                docker : "yirenshao/major_hla_mm_v0.1:latest"  #"liangli92/bmt_v1.5"
                cpu: 32
        }
        output {
                String result = stdout()
                File HLAthena_input_hla = "HLAthena_input_hla.txt"
                File HLAthena_input_peptides = "HLAthena_input_peptides.txt"
                File host_kmers = "host_kmers.txt"
                File donor_kmers = "donor_kmers.txt"

        }
}
### Parse HLA alleles
task parse_alleles_task {
    Int preemptible
    Int memoryGB
    Int diskGB

    File alleles_file
    String alleles_parsed_file =  basename(alleles_file)

    ### using <<< syntax due to parsing problems with {} and awk expression
    command <<<
        # Exit when any command fails
        set -e

        dos2unix ${alleles_file}
        awk '{gsub("HLA-|*","")}1' ${alleles_file} > ${alleles_parsed_file}.parsed.tmp
        awk '{gsub(":","") }1' ${alleles_parsed_file}.parsed.tmp | uniq > ${alleles_parsed_file}.parsed
        rm ${alleles_parsed_file}.parsed.tmp

    >>>

    output {
        File alleles_parsed = '${alleles_parsed_file}.parsed'
    }
    runtime {
        docker: "ssarkizova/hlathena-external:latest"
        memory: "${memoryGB} GB"
        disks: "local-disk ${diskGB} HDD"
        preemptible: 3
    }
}

### Check whether the input is in fasta format and convert to peptides accordingly
task check_for_fasta_and_conver_task {
    Int preemptible
    Int memoryGB
    Int diskGB

    File peptide_list
    String peptide_col_name
    Array[String] lens

    command <<<
        # Exit when any command fails
        set -e

        dos2unix ${peptide_list}
        head -1 ${peptide_list} > ${peptide_list}_header
        if grep -q '>' ${peptide_list}_header;
        then
            echo "FASTA input format detected"
            # Split fasta entries into peptides + context,
            # keep track of the fasta ID each peptide came from
            Rscript /format_input/fasta_reformat.R \
                --infile ${peptide_list} --lens "${sep=' ' lens}" \
                --peptide_col_name ${peptide_col_name} --outfile "peps.txt"
        else
            mv ${peptide_list} "peps.txt"
        fi
    >>>

    output {
        File peptide_list_out = "peps.txt"
    }
    runtime {
        docker: "ssarkizova/hlathena-external:latest"
        memory: "${memoryGB} GB"
        disks: "local-disk ${diskGB} HDD"
        preemptible: "${preemptible}"
    }
}


### Run cleavability predictor
task predict_cleavability_task {
    Int preemptible
    Int memoryGB
    Int diskGB

    File peptide_list
    String peptide_col_name
    String output_file_name = basename(peptide_list)

    command <<<
        # Exit when any command fails
        set -e

        dos2unix ${peptide_list}

        Rscript /clevnn/featurize/featurize_clevnn.R \
            --infile ${peptide_list} \
            --peptide_col_name ${peptide_col_name} \
            --output clevnnfeats.tmp

        python /clevnn/pred/clevnn_pred.py \
            --input_file clevnnfeats.tmp \
            --output_file ${output_file_name}.clev \
            --model_file /clevnn/pred/model_double20_drop0.2_features_375.txt_e10_bs10000_nh150.h5 \
            --features_file /clevnn/pred/features_375.txt

        ### Cleanup
        rm clevnnfeats.tmp
    >>>

    output {
        File peptide_list_clev = "${output_file_name}.clev"
    }
    runtime {
        docker: "ssarkizova/hlathena-external:latest"
        memory: "${memoryGB} GB"
        disks: "local-disk ${diskGB} HDD"
        preemptible: "${preemptible}"
  }
}


### Aggregate per peptide
task aggregate_pep_task {
    Int preemptible
    Int memoryGB
    Int diskGB

    File peptide_list
    String peptide_col_name
    Boolean exists_expr
    String expr_col_name
    String output_file_name = basename(peptide_list)

    command <<<
        # Exit when any command fails
        set -e

        dos2unix ${peptide_list}

        Rscript /aggregate/aggregate_pep.R \
            --peptide_col_name ${peptide_col_name} \
            --exists_expr ${exists_expr} \
            --expr_col_name ${expr_col_name} \
            --infile ${peptide_list} \
            --outfile ${output_file_name}.agg
    >>>

    output {
        File peptide_list_agg = "${output_file_name}.agg"
    }
    runtime {
        docker: "ssarkizova/hlathena-external:latest"
        memory: "${memoryGB} GB"
        disks: "local-disk ${diskGB} HDD"
        preemptible: "${preemptible}"
    }
}


### Ensure expression column is named log2TPM and optionally log-transform
task transform_expression_task {
    Int preemptible
    Int memoryGB
    Int diskGB

    File peptide_list
    String expr_col_name
    Boolean logtransform_expr

    String output_file_name =  basename(peptide_list)

    command <<<
        # Exit when any command fails
        set -e

        dos2unix ${peptide_list}

        Rscript /logtransform_expr/logtransform_expr.R \
            --infile ${peptide_list} \
            --outfile ${output_file_name}.expr \
            --expr_col_name ${expr_col_name} \
            --logtransform_expr "${logtransform_expr}"
    >>>

    output {
        File peptide_list_clev_expr = "${output_file_name}.expr"
    }
    runtime {
        docker: "ssarkizova/hlathena-external:latest"
        memory: "${memoryGB} GB"
        disks: "local-disk ${diskGB} HDD"
        preemptible: "${preemptible}"
    }
}


### Split merged peptides by length
task split_peptides_len_task {
    Int preemptible
    Int memoryGB
    Int diskGB

    File peptide_list
    String peptide_col_name
    Array[String] lens

    command <<<
        # Exit when any command fails
        set -e

        echo peptide_list=${peptide_list}

        dos2unix ${peptide_list}

        for len in ${sep=' ' lens}  ; do
            ### Get the index of the peptide column
            pep_col_id=`head -1 ${peptide_list}  | tr '\t' '\n' | grep -nx ${peptide_col_name} | cut -d: -f1`

            echo pep_col_id=$pep_col_id
            echo len=$len
            echo peps.len.txt=peps.$len.txt

            ### Extract peptides of specified length if any
            if [[ $(tail -n+2 ${peptide_list} | awk -v len="$len" -v pep_col_id="pep_col_id" -F '\t' 'length($'"$pep_col_id"') == len' | wc -l) -ge 1 ]]; then
                head -n 1 ${peptide_list} > peps.$len.txt
                tail -n+2 ${peptide_list} | awk -v len="$len" -v pep_col_id="pep_col_id" -F '\t' 'length($'"$pep_col_id"') == len' >> peps.$len.txt

                echo $len >> lens_present.txt
            fi
        done
    >>>

    output {
        Array[File] pep_len_files = glob('peps.*.txt')
        File lens_present = 'lens_present.txt'
    }
    runtime {
        docker: "ssarkizova/hlathena-external:latest"
        memory: "${memoryGB} GB"
        disks: "local-disk ${diskGB} HDD"
        preemptible: "${preemptible}"
    }
}


### Generate peptide sequences features: dummy
task featurize_encoding_task {
    Int preemptible
    Int memoryGB
    Int diskGB

    File peptide_len_list
    String encoding
    String peptide_col_name

    String output_file_name =  basename(peptide_len_list)

    command <<<
        # Exit when any command fails
        set -e

        echo output_file_name=${output_file_name}
        Rscript /encoding/encoding.R \
            --infile ${peptide_len_list} \
            --outfile ${output_file_name}.${encoding} \
            --peptide_col_name ${peptide_col_name} \
            --coding ${encoding}
    >>>

    output {
        File feature_files = '${output_file_name}.${encoding}'
    }
    runtime {
        docker: "ssarkizova/hlathena-external:latest"
        memory: "${memoryGB} GB"
        disks: "local-disk ${diskGB} HDD"
        preemptible: "${preemptible}"
  }
}


### Predict with MS models (either allele-and-len-specific model, or panpan model when missing)
task predict_ms_allele_task {
    Int preemptible
    Int memoryGB
    Int diskGB
    Int memoryGBx2 = memoryGB*2
    Int diskGBx2 = diskGB*2

    String patient
    File featfile
    Array[String] lens
    Boolean exists_ctex
    Boolean exists_expr
    String allele
    String models_path
    String peptide_col_name
    String features_all
    String feature_sets

    String output_file_name =  basename(featfile)

    command <<<
        # Exit when any command fails
        set -e

        export TMPDIR=/tmp # to fix error: AF_UNIX path too long

        ### Get the len from the input file name
        for infile_len in ${sep=' ' lens} ; do
            if [[ "${featfile}" == *"peps.$infile_len"* ]]; then
                len=$infile_len
            fi
        done
        echo len=$len

        ### Figure out whether to use allele-and-length specific or panpan model
        ### ugly shenanigans splitting grep in two to make wdl happy
        echo grep -w ${allele} /data/model_exists_dat.txt | grep -w $len | cut -d ' ' -f 3
        use_specific=`grep -w ${allele} /data/model_exists_dat.txt | grep -w $len | cut -d ' ' -f 3`
        echo use_specific = $use_specific

        ### Copy models from google bucket and extract the tar.gz file
        if [ "$use_specific" == 1 ]; then
            echo using allele-and-length specific model
            local_model_dir=./${allele} # created by the untar
            gsutil cp ${models_path}${allele}.tar.gz .
            tar -xzf ${allele}.tar.gz .

            models_linear_RDS_file=models_linear_"$len"_${allele}_slim.RDS
            gsutil cp ${models_path}models_linear/$models_linear_RDS_file "$local_model_dir/$models_linear_RDS_file"
        else
            echo using panpan model
            local_model_dir="./panpan"
            mkdir $local_model_dir
            gsutil -m cp -r ${models_path}models_panpan/models_pan_pan_CV/ "$local_model_dir"
            gsutil -m cp -r ${models_path}models_panpan/models_linear_pan_pan/ "$local_model_dir"

            hla=`echo ${allele} | cut -c1` # needed to make WDL happy, $allele:0:1 (with curly) does not work
            models_linear_RDS_file=models_linear_"$len"_"$hla"_slim.RDS # quotes here are very necessary with WDL
            gsutil cp ${models_path}models_panpan/models_linear_pan_pan/"$models_linear_RDS_file" "$local_model_dir/$models_linear_RDS_file"
        fi

        ### Predict MSi
        feature_set=${feature_sets}
        if [ "$use_specific" == 1 ]; then
            ### Execute allele-and-len-specific model
            python /predict/mspredict/ann_pred_py3_optimize.py --sample_name ${patient} --alleles ${allele} --len $len \
                --peptide_col_name ${peptide_col_name} --features_pep_all ${features_all} --feature_set_pep $feature_set \
                --model_path "./" --infile ${featfile} --outpath "./"
        else
            ### Execute panpan model
            python /predict/mspredict/ann_pred_py3_optimize_pan.py --sample_name ${patient} --alleles ${allele} --len $len \
                --peptide_col_name ${peptide_col_name} --features_pep_all ${features_all} --feature_set_pep $feature_set \
                --feature_set_hla_all "allele_feats" --feature_set_hla "allele_feats_PCA" \
                --model_path "./panpan" --infile ${featfile} --outpath "./"
        fi

        ### Predict MSiC, MSiCE, MSiCEB
        if [ "${exists_ctex}" = true ] || [ "${exists_expr}" = true ]; then
            infile=./mspred_${patient}_${allele}_"$len"_${feature_sets}_summary.txt
            outfile=./mspred_${patient}_${allele}_"$len"_${feature_sets}_summary_L.txt
            Rscript /predict/msevallinear/mseval_linear.R \
                    --len $len --allele ${allele} \
                    --models_linear_RDS_file $local_model_dir/$models_linear_RDS_file \
                    --exists_ctex ${exists_ctex} --exists_expr ${exists_expr} \
                    --infile $infile --outfile $outfile
            mv $outfile $infile
        fi

        ### Cleanup
        if [ -f ./${allele}.tar.gz ]; then
            rm ${allele}.tar.gz
            rm -r ./${allele}
        else
            rm -r ./panpan
        fi
    >>>

    output {
        Array[File] mspred_files = glob('*_summary.txt')
    }
    runtime {
        docker: "ssarkizova/hlathena-external:latest"
        memory: "${memoryGBx2} GB"
        disks: "local-disk ${diskGBx2} HDD"
        preemptible: "${preemptible}"
  }
}


### Merge MS predictions
task merge_mspreds_task {
    Int preemptible
    Int memoryGB
    Int diskGB
    Int memoryGBx2 = memoryGB*2
    Int diskGBx2 = diskGB*2

    String patient
    Array[String] alleles
    String len
    Boolean exists_ctex
    Boolean exists_expr
    String peptide_col_name
    Array[Array[File]] mspred_files

    command <<<
        ### Get the index of the peptide column
        echo peptide_col_name=${peptide_col_name}
        echo len=${len}
        output_file_name="mspreds.${patient}.${len}"
        echo output_file_name=$output_file_name

        ### Split the Array[Array[File]] of input files - There should be a better way to do this
        mspred_files_str=""
        for mspred_files_inner in ${sep=' ' mspred_files} ; do
            echo mspred_files_inner=$mspred_files_inner
            for mspred_file in $(echo $mspred_files_inner | tr "[,]" "\n") ; do
                echo mspred_file=$mspred_file
                if [ ! -z $mspred_file ] ; then
                    if [ -z $mspred_files_str ] ; then
                        mspred_files_str=$mspred_file
                    else
                        mspred_files_str=$mspred_files_str,$mspred_file
                    fi
                fi
            done
        done
        echo $mspred_files_str

        Rscript /predict/msmerge/msmerge_optimized.R --patient ${patient} --alleles ${sep=',' alleles} --len ${len} \
                    --exists_ctex ${exists_ctex} --exists_expr ${exists_expr} --peptide_col_name ${peptide_col_name} \
                    --infiles $mspred_files_str --outfile $output_file_name \
                    --run_netmhc "false" --netmhc_EL "" --netmhc_BA ""
    >>>

    output {
        File mspreds_file_long = "mspreds.${patient}.${len}.long.txt"
        File mspreds_file_wide = "mspreds.${patient}.${len}.wide.txt"
    }
    runtime {
        docker: "ssarkizova/hlathena-external:latest"
        memory: "${memoryGBx2} GB"
        disks: "local-disk ${diskGBx2} HDD"
        preemptible: "${preemptible}"
  }
}


task get_ranks_task {
    Int preemptible
    Int memoryGB
    Int diskGB
    Int memoryGBx2 = memoryGB*2
    Int diskGBx2 = diskGB*2

    String patient
    Array[String] alleles
    Array[String] lens
    Boolean exists_ctex
    Boolean exists_expr
    String models_path
    String peptide_col_name
    File preds_file
    String preds_file_name = basename(preds_file)

    command <<<
        # Exit when any command fails
        set -e

        export TMPDIR=/tmp # to fix error: AF_UNIX path too long

        ### Get the len from the input file name
        for infile_len in ${sep=' ' lens} ; do
            if [[ "${preds_file}" == *"$infile_len.wide.txt"* ]]; then
                len=$infile_len
            fi
        done
        echo len=$len

        ### Copy ecdf files from google bucket
        for allele in ${sep=' ' alleles}; do
            ### Figure out whether to use allele-and-length specific or panpan model
            ### ugly shenanigans splitting grep in two to make wdl happy
            echo grep -w "$allele" /data/model_exists_dat.txt | grep -w $len | cut -d ' ' -f 3
            use_specific=`grep -w "$allele" /data/model_exists_dat.txt | grep -w $len | cut -d ' ' -f 3`
            echo use_specific = $use_specific

            if [ "$use_specific" == 1 ]; then
                echo "$allele": using allele-and-length specific model
                model_ecdf=ecdf_panpan_"$allele"_nmhc_mhcflurry.RDS # 'panpan' here is just a typo, these are specific ecdfs
                gsutil cp ${models_path}ecdf/"$model_ecdf" .
            else
                echo "$allele": using panpan model
                model_ecdf=ecdf_panpan_"$allele".RDS
                gsutil cp ${models_path}models_panpan/ecdf/"$model_ecdf" .
            fi
        done

        ### Ranks
        Rscript /predict/ecdf/ecdf_wide.R \
                    --len $len --alleles ${sep=',' alleles} \
                    --exists_ctex ${exists_ctex} --exists_expr ${exists_expr} --run_netmhc "false" \
                    --ecdf_model_path "./" \
                    --infile ${preds_file} --outfile ${preds_file_name}.ranks
    >>>

    output {
        File ranks_file = "${preds_file_name}.ranks"
    }
    runtime {
        docker: "ssarkizova/hlathena-external:latest"
        memory: "${memoryGBx2} GB"
        disks: "local-disk ${diskGBx2} HDD"
        preemptible: "${preemptible}"
  }
}


### Concatenate ranks files for all lengths, assign alleles, provide summaries
task ranks_concat_lens_task {
    Int preemptible
    Int memoryGB
    Int diskGB
    Int diskGBx2 = diskGB*2

    String patient
    String peptide_col_name
    Boolean exists_ctex
    Boolean exists_expr
    String assign_by_ranks_or_scores
    String assign_threshold
    String assign_colors
    Array[File] ranks_len_files
    String outfile = "${patient}.ranks.cat"

    ### using <<< syntax due to parsing problems with {} and awk expression
    command <<<
        # Exit when any command fails
        set -e

        # Only keep header from first file:
        #    FNR is the number of lines (records) read so far in the current file.
        #     NR is the number of lines read overall.
        awk 'FNR==1 && NR!=1{next;}{print}' ${sep=" " ranks_len_files} > ${outfile}

        # Assign peptides to alleles and make summary plots
        Rscript /predict/assign/assign.R \
            --sample_name ${patient} --peptide_col_name ${peptide_col_name} \
            --exists_ctex ${exists_ctex} --exists_expr ${exists_expr} --run_netmhc "false" \
            --assign_by_ranks_or_scores ${assign_by_ranks_or_scores} --assign_threshold ${assign_threshold} --assign_colors "${assign_colors}" \
            --infile ${outfile} --outpath "./"
    >>>

    output {
        File sample_predictions = "${patient}_predictions.txt"
        File sample_allele_assignment_counts = "${patient}_allele_assignment_counts.txt"
        File sample_allele_assignment_counts_for_plot = "${patient}_allele_assignment_counts_for_plot.txt"
        File sample_plots = "${patient}.pdf"
    }
    runtime {
        docker: "ssarkizova/hlathena-external:latest"
        memory: "${memoryGB} GB"
        disks: "local-disk ${diskGBx2} HDD"
        preemptible: "${preemptible}"
  }
}