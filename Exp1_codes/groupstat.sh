## transform hires to the MNI template
subs="01 03 05 07 09 11 13 15 17 19 21 23"
for subID in $subs
do
   cd ${subID}/Img_data/MNI_transform
   
   ## copy
   cp ../T1_al2Surf+orig.* ./
   cp ../${subID}Avg_fullCross_MSTaskLoc.results/stats.${subID}Avg_MSTaskLoc2+orig.* ./
   cp ../${subID}Avg_fullCross_MSTaskLoc.results/stats.${subID}Avg_fullCross_MSTaskLoc+orig.* ./
   
   # @auto_tlrc -base TT_avg152T1+tlrc -suffix _mni -input T1_al2Surf+orig
   @auto_tlrc -base MNI_avg152T1+tlrc  -suffix _mni -init_xform AUTO_CENTER -input T1_al2Surf+orig. -dxyz 1
   @auto_tlrc -apar T1_al2Surf_mni+tlrc -input stats.${subID}Avg_MSTaskLoc2+orig -suffix _mni -dxyz 3
   @auto_tlrc -apar T1_al2Surf_mni+tlrc -input stats.${subID}Avg_fullCross_MSTaskLoc+orig -suffix _mni -dxyz 3
 
done

## t-test 
3dttest++ -setA 'stats.01Avg_MSTaskLoc2_mni+tlrc.HEAD[1]' 'stats.03Avg_MSTaskLoc2_mni+tlrc.HEAD[1]' \
      'stats.05Avg_MSTaskLoc2_mni+tlrc.HEAD[1]' 'stats.07Avg_MSTaskLoc2_mni+tlrc.HEAD[1]' \
      'stats.09Avg_MSTaskLoc2_mni+tlrc.HEAD[1]' 'stats.11Avg_MSTaskLoc2_mni+tlrc.HEAD[1]' \
      'stats.13Avg_MSTaskLoc2_mni+tlrc.HEAD[1]' 'stats.15Avg_MSTaskLoc2_mni+tlrc.HEAD[1]' \
      'stats.17Avg_MSTaskLoc2_mni+tlrc.HEAD[1]' 'stats.19Avg_MSTaskLoc2_mni+tlrc.HEAD[1]' \
      'stats.21Avg_MSTaskLoc2_mni+tlrc.HEAD[1]' 'stats.23Avg_MSTaskLoc2_mni+tlrc.HEAD[1]' \
-setB 'stats.01Avg_MSTaskLoc2_mni+tlrc.HEAD[4]' 'stats.03Avg_MSTaskLoc2_mni+tlrc.HEAD[4]' \
      'stats.05Avg_MSTaskLoc2_mni+tlrc.HEAD[4]' 'stats.07Avg_MSTaskLoc2_mni+tlrc.HEAD[4]' \
      'stats.09Avg_MSTaskLoc2_mni+tlrc.HEAD[4]' 'stats.11Avg_MSTaskLoc2_mni+tlrc.HEAD[4]' \
      'stats.13Avg_MSTaskLoc2_mni+tlrc.HEAD[4]' 'stats.15Avg_MSTaskLoc2_mni+tlrc.HEAD[4]' \
      'stats.17Avg_MSTaskLoc2_mni+tlrc.HEAD[4]' 'stats.19Avg_MSTaskLoc2_mni+tlrc.HEAD[4]' \
      'stats.21Avg_MSTaskLoc2_mni+tlrc.HEAD[4]' 'stats.23Avg_MSTaskLoc2_mni+tlrc.HEAD[4]' \
-prefix TaskOn_Off -paired

## average all subs' hi-res.
cd group_stats

3dcalc -a 01_T1_mni+tlrc -prefix AvgfullCross_T1_mni -b 03_T1_mni+tlrc -c 05_T1_mni+tlrc -d 07_T1_mni+tlrc -e 09_T1_mni+tlrc \
-f 11_T1_mni+tlrc -g 13_T1_mni+tlrc -h 15_T1_mni+tlrc -i 17_T1_mni+tlrc -j 19_T1_mni+tlrc -k 21_T1_mni+tlrc -l 23_T1_mni+tlrc \
-exp 'mean(a,b,c,d,e,f,g,h,i,j,k,l)'

