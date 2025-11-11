# Attempt of a script to perform vertex-wise analyses on two neuroimage metrics 

# Usage: 3 elements needed:
#   - map2map folder with FSGD and Contrasts subfolders
#   - FSGD files (map2map_FSGD_Var.txt, in FSGD folder, with single class map2map)
#   - Contrast files (contrasts, in Contrasts folder)

export SUBJECTS_DIR=/home/oriol/Desktop/TFM/donsurf/donsurf_TP1/derivatives
WDIR=/home/oriol/Desktop/TFM/donsurf/3rd_analyses/map2map

METRIC=("thickness GM-MD-donsurf1.0.0")
hemisphere=("lh rh")
Smoothing=("15")
contrasts=("contrasts")
Variables=("map2map_no_covar")

echo %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Preparing files and folders %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

############# Folder structure 
for Var in $Variables; do
    mkdir $WDIR/$Var/sig.cluster_screenshot
    for cxt in $contrasts; do
        mkdir $WDIR/$Var/sig.cluster_screenshot/$cxt;
    done;
    mkdir $WDIR/$Var/concat;
    mkdir $WDIR/$Var/glmfit;
    mkdir $WDIR/$Var/tmp;
done

############# FSGD to unix to .mat
for Var in $Variables; do
    cat $WDIR/$Var/FSGD/map2map_FSGD_${hemi}.fsgd | sed 's/\r/\n/g' > $WDIR/$Var/FSGD/map2map_FSGD_${hemi}2.fsgd
    tr ',' '.' < $WDIR/$Var/FSGD/map2map_FSGD_${Var}.txt > $WDIR/$Var/tmp/map2map_FSGD_points.txt
    tr '\t' ' ' < $WDIR/$Var/tmp/map2map_FSGD_points.txt > $WDIR/$Var/FSGD/map2map_FSGD_${Var}.fsgd
    dos2unix $WDIR/$Var/FSGD/map2map_FSGD_${Var}.fsgd
    sed 's/map2map/0/g; /Class\|Title\|Descriptor\|Variable/d; s/Input[^ ]* *\|sub[^ ]* *//g' $WDIR/$Var/FSGD/map2map_FSGD_${Var}.fsgd > $WDIR/$Var/tmp/pre_X.mat
    sed 's/0 /1 /g' $WDIR/$Var/tmp/pre_X.mat > $WDIR/$Var/glmfit/X.mat
    wc -l $WDIR/$Var/glmfit/X.mat > $WDIR/$Var/tmp/number
    awk '{print $1}' $WDIR/$Var/tmp/number > $WDIR/$Var/tmp/row_number
    for i in $(seq 1 `cat $WDIR/$Var/tmp/row_number`); 
           do echo 1; 
    done > $WDIR/$Var/glmfit/X2.mat  # matrix without covars
done

########### Subjects lists

for Var in $Variables; do
    grep -o "sub-......" $WDIR/$Var/FSGD/map2map_FSGD_${Var}.fsgd > $WDIR/$Var/tmp/pre_sub_list.txt
    sed 's/ .*//g' $WDIR/$Var/tmp/pre_sub_list.txt > $WDIR/$Var/concat/sub_list.txt
    for metr in $METRIC; do
        for hemi in $hemisphere; do
            for smth in $Smoothing; do
              for sub in `cat $WDIR/$Var/concat/sub_list.txt`; 
              do echo $SUBJECTS_DIR/$sub/surf/${hemi}.${metr}.fwhm${smth}.fsaverage.mgh; 
              done > $WDIR/$Var/concat/${hemi}.${metr}.${smth}.path_list.txt;
            done
        done
    done
    rm $WDIR/$Var/tmp -r;
done



echo %%%%%%%%%%%%%%%%%%%%%%%%%%%%% Concatenating images and creating masks %%%%%%%%%%%%%%%%%%%%%%%%%%

########## Concatenating images: mri_concat

for Var in $Variables; do
    for metr in $METRIC; do
                for hemi in $hemisphere; do
                        for smth in $Smoothing; do
                        echo Concatenating Metric=$metr; echo Hemisphere=$hemi; echo Smoothing=$smth;
                        mri_concat --f $WDIR/$Var/concat/${hemi}.${metr}.${smth}.path_list.txt \
                        --o $WDIR/$Var/concat/${hemi}_${metr}_${smth}_concat.mgh;
                        done
                done
    done
done

########## concenating images: mris_preproc (alternatively)

# for metr in $METRIC; do
#            for hemi in $hemisphere; do
#                    for smth in $Smoothing; do
#                    echo concatenating Metric=$metr; echo Hemisphere=$hemi; echo Smoothing=$smth;
#                    mris_preproc --fsgd $WDIR/$Var/FSGD/map2map_FSGD_${Var}.fsgd \
#                    --cache-in ${metr}.fwhm${smth}.fsaverage \
#                    --target fsaverage --hemi $hemi \
#                    --out $WDIR/$Var/concat/${hemi}_${metr}_${smth}_concat.mgh
#                    done
#            done
# done

########## Generating masks (only for cMD)
  
for Var in $Variables; do
    for hemi in $hemisphere; do
        for smth in $Smoothing; do
            echo Hemisphere=$hemi; echo Smoothing=$smth;
            mri_binarize --i $WDIR/$Var/concat/${hemi}_GM-MD-donsurf1.0.0_${smth}_concat.mgh \
            --min .0000001 \
            --o $WDIR/$Var/concat/${hemi}.${smth}.framemask.mgh;
        done
    done
done

echo %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Running MRI_GLMFIT %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

########### mri_glmfit 

for Var in $Variables; do
    for cxt in $contrasts; do
        for hemi in $hemisphere; do
            for smth in $Smoothing; do
                echo Hemisphere=$hemi; echo Smoothing=$smth;
                mri_glmfit --glmdir $WDIR/$Var/glmfit/${hemi}.${smth}.output_glmdir \
                           --y $WDIR/$Var/concat/${hemi}_thickness_${smth}_concat.mgh \
                           --fsgd $WDIR/$Var/FSGD/map2map_FSGD_${Var}.fsgd \
                           --C $WDIR/$Var/Contrasts/$cxt \
                           --surf fsaverage $hemi \
                           --frame-mask $WDIR/$Var/concat/${hemi}.${smth}.framemask.mgh \
                           --cortex \
                           --no-prune \
                           --pvr $WDIR/$Var/concat/${hemi}_GM-MD-donsurf1.0.0_${smth}_concat.mgh;    
            done
        done
    done
done
        
############ Cluster correction with mri_glmfit-sim

export LC_NUMERIC="en_US.UTF-8"

for Var in $Variables; do
    for hemi in $hemisphere; do
          for smth in $Smoothing; do
               echo Hemisphere=$hemi; echo Smoothing=$smth;
               mri_glmfit-sim --glmdir $WDIR/$Var/glmfit/${hemi}.${smth}.output_glmdir --cache 1.3 abs --cwp 0.05 --2spaces
          done
    done
done

#echo %%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Capturing Freeview Images %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

############ freeview screenshots

#for Var in $Variables; do
#    for cxt in $contrasts; do
#        for hemi in $hemisphere; do
#            for smth in $Smoothing; do
#                echo Variable=$Var; echo Hemisphere=$hemi; echo Smoothing=$smth;
#                freeview -f $SUBJECTS_DIR/fsaverage/surf/${hemi}.pial_semi_inflated:overlay=$WDIR/$Var/glmfit/${hemi}.${smth}.output_glmdir/${cxt}/cache.th13.abs.sig.cluster.mgh:overlay_threshold=2.3,3:overlay_frame=1 -viewport 3d -layout 1 -nocursor -colorscale --ss $WDIR/$Var/sig.cluster_screenshot/$cxt/${hemi}.${smth}_left.png 
#                            freeview -f $SUBJECTS_DIR/fsaverage/surf/${hemi}.pial_semi_inflated:overlay=$WDIR/$Var/glmfit/${hemi}.${smth}.output_glmdir/${cxt}/cache.th13.abs.sig.cluster.mgh:overlay_threshold=2.3,3:overlay_frame=1 -viewport 3d -layout 1 -nocursor -colorscale -cam azimuth 180 --ss $WDIR/$Var/sig.cluster_screenshot/$cxt/${hemi}.${smth}_right.png   
#            done
#        done
#    done
#done
echo Done!
