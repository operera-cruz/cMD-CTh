#!/usr/bin/env bash

# Attempt of a script to plot effect sizes on significant clusters

# Usage: specify SUBJECTS_DIR (subject data) and WDIR (analysis folder). Run with --help for other options.

export SUBJECTS_DIR=/home/oriol/Desktop/TFM/donsurf/donsurf_TP1/derivatives
WDIR=/home/oriol/Desktop/TFM/donsurf/3rd_analyses/map2map

METRIC=("GM-MD-donsurf1.0.0 thickness")
group=("PACC_AG1 PACC_AG2 PACC_AG3")
hemisphere=("lh rh")
Smoothing=("15")
Variables=("map2map_no_covar")

# Parse flags and warnings

CLEAN=false
clus="surfcluster"
FV=false
thmin=""
thmax=""

while [[ $# -gt 0 ]]; do
  case "$1" in
        --clean)
            CLEAN=true
            ;;
        --th-min)
            shift
            thmin="$1"
            ;;
        --th-max)
            shift
            thmax="$1"
            ;;
        --clus=*)
            clus="${1#*=}"
            ;;
        --fv)
            FV=true
            ;;
        --help)
            echo "Usage: $0 
            Must specify within the script: SUBJECTS_DIR, WDIR, METRIC, hemisphere, Smoothing, Variables
            [--clean] will remove tmp files
            [--raw-pval] will extract masks from raw sig.mgh files assuming pval>.05
            [--fv] will capture freeview screenshots
            [--th-min] [--th-max] to define overlay max threshold in freeview screenshots 
            [--clus to specify output sig file to extract effect-size from. Options are:
                surfcluster --> for mri_surfcluster sig.cluster.mgh files (default)
                rawpval     --> for mri_surfcluster sig.mgh files (pval<.05)
                glmfitsim   --> for glmfit-sim sig.cluster.mgh files  "
            exit 0
            ;;
        *)
            echo "Unknown option: $arg"
            sleep 2
            ;;
    esac
    shift
done

if [ "$clus" = "surfcluster" ]; then
    echo "Running effect-size calculation for mri_surfcluster output (default) with masks based on cluster-corrected stats"
    sleep 2
elif [ "$clus" = "rawpval" ];then
    echo "Running effect-size calculation for mri_surfcluster output but BEWARE! Masks obtained from raw stats (pval< .05)"
    sleep 2
elif [ "$clus" = "glmfitsim" ];then
    echo "Running effect-size calculations for glmfit-sim output"
    sleep 2
else
    echo "Error: Unkown --clus option 
    options: surfcluster, rawpval, glmfitsim
    Exiting..."
    sleep 2
    exit 1
fi

# Set folders
for met in $METRIC; do
#for grp in $group; do
    for Var in $Variables; do
        mkdir $WDIR/$Var/effect_size;
        mkdir $WDIR/$Var/effect_size/tmp
    done
#done
done

######### Effect-size calculation Steps: 
#    1) extract frames (only for mri_surfcluster output)
#    2) binarize clusters 
#        for negative and positive values sepparately
#        then unify them
#    3) multiply mask by effect-size (gamma=C*B)

for met in $METRIC; do
#for grp in $group; do
    for Var in $Variables; do
        for hemi in $hemisphere; do
            for smth in $Smoothing; do
            if [ "$clus" = "surfcluster" ]||[ "$clus" = "rawpval" ]; then
                echo "%%%%%%%%%%%%%%%%%%%%%%% Extracting frames %%%%%%%%%%%%%%%%%%%%%%%"
                mri_convert $WDIR/$Var/GLMDIR/${hemi}_results/contrasts/mc-z.abs.th13.sig.cluster.mgh --frame 2 $WDIR/$Var/effect_size/tmp/${hemi}.mc-z.abs.th13.sig.cluster.f1.mgh
                echo "%%%%%%%%%%%%%%%%%%%%%%% Generating masks %%%%%%%%%%%%%%%%%%%%%%%";
                    if [ "$clus" = "surfcluster" ]; then
                        mri_binarize --i $WDIR/$Var/effect_size/tmp/${hemi}.mc-z.abs.th13.sig.cluster.f1.mgh \
                        --min -inf --max -1.30103 --o $WDIR/$Var/effect_size/tmp/${hemi}.${smth}.neg_mask.mgh;
                        mri_binarize --i $WDIR/$Var/effect_size/tmp/${hemi}.mc-z.abs.th13.sig.cluster.f1.mgh \
                        --min 1.30103 --o $WDIR/$Var/effect_size/tmp/${hemi}.${smth}.pos_mask.mgh;
                        mris_calc -o $WDIR/$Var/effect_size/tmp/${hemi}.${smth}.full_mask.mgh $WDIR/$Var/effect_size/tmp/${hemi}.${smth}.neg_mask.mgh add $WDIR/$Var/effect_size/tmp/${hemi}.${smth}.pos_mask.mgh;
                     else
                        mri_binarize --i $WDIR/$Var/GLMDIR/${hemi}_results/contrasts/sig.mgh \
                        --min -inf --max -1.30103 --o $WDIR/$Var/effect_size/tmp/${hemi}.${smth}.neg_mask.ncc.mgh;
                        mri_binarize --i $WDIR/$Var/GLMDIR/${hemi}_results/contrasts/sig.mgh \
                        --min 1.30103 --o $WDIR/$Var/effect_size/tmp/${hemi}.${smth}.pos_mask.ncc.mgh;
                        mris_calc -o $WDIR/$Var/effect_size/tmp/${hemi}.${smth}.full_mask.mgh $WDIR/$Var/effect_size/tmp/${hemi}.${smth}.neg_mask.ncc.mgh add $WDIR/$Var/effect_size/tmp/${hemi}.${smth}.pos_mask.ncc.mgh;
                     fi
                 echo " %%%%%%%%%%%%%%%%%%%%%%% Generating effect-size maps for $Var $hemi with $smth smoothing %%%%%%%%%%%%%%%%%%%%%%%"
                        mris_calc -o $WDIR/$Var/effect_size/${hemi}.${smth}.effect_size_map.mgh $WDIR/$Var/effect_size/tmp/${hemi}.${smth}.full_mask.mgh mul $WDIR/$Var/GLMDIR/${hemi}_results/contrasts/gamma.mgh
            elif [ "$clus" = "glmfitsim" ];then
                echo "%%%%%%%%%%%%%%%%%%%%%%% Generating masks %%%%%%%%%%%%%%%%%%%%%%%";
                        mri_binarize --i $WDIR/$Var/glmfit/${hemi}.${smth}.output_glmdir/contrasts/cache.th13.abs.sig.cluster.mgh \
                        --min -inf --max -1.30103 --o $WDIR/$Var/effect_size/tmp/${hemi}.${smth}.neg_mask.sim.mgh;
                        mri_binarize --i $WDIR/$Var/glmfit/${hemi}.${smth}.output_glmdir/contrasts/cache.th13.abs.sig.cluster.mgh \
                        --min 1.30103 --o $WDIR/$Var/effect_size/tmp/${hemi}.${smth}.pos_mask.sim.mgh;
                        mris_calc -o $WDIR/$Var/effect_size/tmp/${hemi}.${smth}.full_mask.mgh $WDIR/$Var/effect_size/tmp/${hemi}.${smth}.neg_mask.sim.mgh add $WDIR/$Var/effect_size/tmp/${hemi}.${smth}.pos_mask.sim.mgh;
                echo "%%%%%%%%%%%%%%%%%%%%%%% Generating effect-size maps for $Var $hemi with $smth smoothing %%%%%%%%%%%%%%%%%%%%%%%"
                        mris_calc -o $WDIR/$Var/effect_size/${hemi}.${smth}.effect_size_map.mgh $WDIR/$Var/effect_size/tmp/${hemi}.${smth}.full_mask.mgh mul $WDIR/$Var/glmfit/${hemi}.${smth}.output_glmdir/contrasts/gamma.mgh
            fi
                    
                    
### freeview screenshots
                    if $FV; then
                        if [ -z "$thmin" ] && [ -z "$thmax" ]; then
                            echo " %%%%%%%%%%%%%%%%%%%%%%% Capturing images for $hemi $met x $Var with default overaly threshold (.A) %%%%%%%%%%%%%%%%%%%%%%%"
                            freeview -f $SUBJECTS_DIR/fsaverage/surf/${hemi}.pial_semi_inflated:overlay=$WDIR/$Var/effect_size/${hemi}.${smth}.effect_size_map.mgh:overlay_frame=1 -viewport 3d -layout 1 -nocursor -colorscale --ss $WDIR/$Var/effect_size/${hemi}.${smth}_lateral.A.png 
                            freeview -f $SUBJECTS_DIR/fsaverage/surf/${hemi}.pial_semi_inflated:overlay=$WDIR/$Var/effect_size/${hemi}.${smth}.effect_size_map.mgh:overlay_frame=1 -viewport 3d -layout 1 -nocursor -colorscale -cam azimuth 180 --ss $WDIR/$Var/effect_size/${hemi}.${smth}_medial.A.png  
                        elif [ -n "$thmin" ] && [ -n "$thmax" ]; then
                            echo "%%%%%%%%%%%%%%%%%%%%%%% Capturing images for $hemi $met x $Var with overlay threshold set to min $thmin max $thmax (.B) %%%%%%%%%%%%%%%%%%%%%%%"
                            freeview -f $SUBJECTS_DIR/fsaverage/surf/${hemi}.pial_semi_inflated:overlay=$WDIR/$Var/effect_size/${hemi}.${smth}.effect_size_map.mgh:overlay_threshold="$thmin","$thmax","$TH":overlay_frame=1 -viewport 3d -layout 1 -nocursor -colorscale --ss $WDIR/$Var/effect_size/${hemi}.${smth}_lateral.B.png 
                            freeview -f $SUBJECTS_DIR/fsaverage/surf/${hemi}.pial_semi_inflated:overlay=$WDIR/$Var/effect_size/${hemi}.${smth}.effect_size_map.mgh:overlay_threshold="$thmin","$thmax":overlay_frame=1 -viewport 3d -layout 1 -nocursor -colorscale -cam azimuth 180 --ss $WDIR/$Var/effect_size/${hemi}.${smth}_medial.B.png
                        else 
                            echo "Error: both --th-min and --th-max must be provided together" >&2
                            exit 1
                        fi
                    else
                        echo "Done! To Capture freeview screenshots use --fv (optional --th-min --th-max)"
                    fi
            done
        done
    done
#done
done

###### remove tmp files if --clean
for met in $METRIC; do
    for Var in $Variables; do
        if [[ "$CLEAN" == true ]]; then
        echo "  Removing intermediate positive and negative masks"
            rm $WDIR/$Var/effect_size/tmp
        fi
    done
done


