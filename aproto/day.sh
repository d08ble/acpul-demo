#!/bin/sh


convert -font Ravie -pointsize 72  label:'1' -border 10 \
          -tile tile_aqua.jpg   -draw "color 0,0 reset"  \
          -tile tile_water.jpg -gravity center -annotate +0+0 'Get Wet!' \
          autosize_wet.jpg