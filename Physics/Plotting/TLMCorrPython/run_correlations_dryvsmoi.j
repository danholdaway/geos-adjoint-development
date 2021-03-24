# Run the correlations script for all the variables


./correlations_dryvsmoist_GQ1.py --end_date 2013012800 --root_tmpl $ARCHIVE -s --fhr 24 --region tropics --zaxis linear --variable sphu
mv 2013010100-2013012800-f24_tropics_PH2_sphu.pdf 2013010100-2013012800-f24_tropics_dryvsmoi_GQ1_sphu.pdf

./correlations_dryvsmoist_GQ1.py --end_date 2013012800 --root_tmpl $ARCHIVE -s --fhr 48 --region tropics --zaxis linear --variable sphu
mv 2013010100-2013012800-f48_tropics_PH2_sphu.pdf 2013010100-2013012800-f48_tropics_dryvsmoi_GQ1_sphu.pdf

./correlations_dryvsmoist_GQ2.py --end_date 2013012800 --root_tmpl $ARCHIVE -s --fhr 24 --region tropics --zaxis linear --variable sphu
mv 2013010100-2013012800-f24_tropics_PH2_sphu.pdf 2013010100-2013012800-f24_tropics_dryvsmoi_GQ2_sphu.pdf

./correlations_dryvsmoist_GQ2.py --end_date 2013012800 --root_tmpl $ARCHIVE -s --fhr 48 --region tropics --zaxis linear --variable sphu
mv 2013010100-2013012800-f48_tropics_PH2_sphu.pdf 2013010100-2013012800-f48_tropics_dryvsmoi_GQ2_sphu.pdf


rm *.eps
rm *.png



