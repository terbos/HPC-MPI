#!/bin/sh
cd ./OUT
rm wave*.ps
for i in *
do
	cat ../head.ps $i ../tail.ps > wave"$i".ps
done
convert -delay 20 wave*.ps ../animation.gif
