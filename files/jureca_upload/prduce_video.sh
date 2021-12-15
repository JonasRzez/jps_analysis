python ini_mc_large_i.py
python mc_l_1.py
#wait

#cd exp_results/plots/video_test/

#rm *

#cd ../../..

python disc_plot.py

#cd exp_results/plots/video_test/



#ffmpeg -f image2 -r 5 -i %03d.png -vf "pad=ceil(iw/2)*2:ceil(ih/2)*2" -pix_fmt yuv420p  i23_5fps.mp4

#open i23_5fps.mp4
