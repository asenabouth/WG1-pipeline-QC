#!/usr/local/envs/py36/bin python3


rule gtm_preprocess:
    input:
        mhc = input_dict["pipeline_dir"] + "/MHC_location.txt",
        bed = output_dict["output_dir"] + "/hrc_check/strand_check-updated-chr{chr}.bed"
    output:
        output_dict["output_dir"] + "/gtm_preprocess/gtm_preprocess_chr{chr}_prune.prune.in"
    resources:
        mem_per_thread_gb=lambda wildcards, attempt: attempt * GTM_check_dict["gtm_preprocess_memory"],
        disk_per_thread_gb=lambda wildcards, attempt: attempt * GTM_check_dict["gtm_preprocess_memory"]
    threads: GTM_check_dict["gtm_preprocess_threads"]
    params:
        bed = output_dict["output_dir"] + "/hrc_check/strand_check-updated-chr{chr}",
        out = output_dict["output_dir"] + "/gtm_preprocess/gtm_preprocess_chr{chr}",
        hrc_check = output_dict["output_dir"] + "/hrc_check",
        sif = input_dict["singularity_image"],
        out_base = output_dict["output_dir"] + "/gtm_preprocess/gtm_preprocess"
    shell:
        """
        if [[ {wildcards.chr} == 6 ]]
        then
            singularity exec {params.sif} plink --exclude {input.mhc} --bfile {params.bed} --indep-pairwise 1000 10 0.02 --maf 0.05 --out {params.out} --make-bed
        else
        then
            singularity exec {params.sif} plink --bfile {params.bed} --indep-pairwise 1000 10 0.02 --maf 0.05 --out {params.out} --make-bed
        fi
        """

rule gtm_prune:
    input:
        bed = output_dict["output_dir"] + "/hrc_check/strand_check-updated-chr{chr}.bed",
        prune = output_dict["output_dir"] + "/gtm_preprocess/gtm_preprocess_chr{chr}.prune.in"
    output:
        bed = output_dict["output_dir"] + "/gtm_prune/gtm_prune_chr{chr}.bed"
    resources:
        mem_per_thread_gb = lambda wildcards, attempt: attempt * GTM_check_dict["gtm_prune_memory"],
        disk_per_thread_gb = lambda wildcards, attempt: attempt * GTM_check_dict["gtm_prune_memory"]
    threads: GTM_check_dict["gtm_prune_threads"]
    params:
        bed = output_dict["output_dir"] + "/hrc_check/strand_check-updated-chr{chr}",
        sif = input_dict["singularity_image"],
        out = output_dict["output_dir"] + "/gtm_prune/gtm_prune_chr{chr}"
    shell:
        """
        singularity exec {params.sif} plink --bfile {params.bed} --extract {input.prune} --out {params.out} --make-bed
        """

rule gtm_merge:
    input:
        output_dict["output_dir"] + "/gtm_prune/gtm_prune_chr{chr}.bed"
    output:
        merge = output_dict["output_dir"] + "gtm_merge/tomerge.txt",
        out = output_dict["output_dir"] + "/gtm_merge/gtm_merge_chr.raw",
        mat = output_dict["output_dir"] + "/gtm_merge/gtm_merge_chr.mat"
    resources:
        mem_per_thread_gb=lambda wildcards, attempt: attempt * GTM_check_dict["gtm_merge_memory"],
        disk_per_thread_gb=lambda wildcards, attempt: attempt * GTM_check_dict["gtm_merge_memory"]
    threads: GTM_check_dict["gtm_merge_threads"]
    params:
        sif = input_dict["singularity_image"],
        out = output_dict["output_dir"] + "/gtm_merge/gtm_merge_chr"
    shell:
        """
        singularity exec {params.sif} echo {input} | singularity exec {params.sif} tr '.bed' '\n' > {output.merge}
        singularity exec {params.sif} plink --merge-list {output.merge} --recodeA --out {params.out}
        singularity exec {params.sif} tail -n +2 {outut.out} > {output.mat}
        """

rule gtm_projection:
    input:
        mat = output_dict["output_dir"] + "/gtm_merge/gtm_merge_chr.mat"
    output:
        output_dict["output_dir"] + "/gtm_projection/gtm_projection_1000G_indiv_predictions.csv",
        output_dict["output_dir"] + "/gtm_projection/gtm_projection_1000G_indiv_probabilities.csv",
        output_dict["output_dir"] + "/gtm_projection/gtm_projection_1000G_group_probabilities.csv"
    resources:
        mem_per_thread_gb=lambda wildcards, attempt: attempt * GTM_check_dict["gtm_projection_memory"],
        disk_per_thread_gb=lambda wildcards, attempt: attempt * GTM_check_dict["gtm_projection_memory"]
    threads: GTM_check_dict["gtm_projection_threads"]
    params:
        sif = input_dict["singularity_image"],
        out_base = output_dict["output_dir"] + "/gtm_projection/gtm_projection"
    shell:
        """
        singularity exec {params.sif} python runGTM.py \
            --model GTM \
            --data /opt/ancestry_viz/recoded_1000G.noadmixed.mat \
            --test {input.mat} \
            --labels /opt/ancestry_viz/recoded_1000G.raw.noadmixed.lbls3_3 \
            --labeltype discrete \
            --out gtm_projection_1000G \
            --pca \
            --n_components 10 \
            --missing \
            --missing_strategy median \
            --random_state 8 
        """