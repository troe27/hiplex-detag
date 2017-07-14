#!/bin/bash
python ../detag_v3.py --tags ./tags_no_CAA.txt --heel ./heel.txt --pairedfastq --forward first_100k_R1.fastq --reverse first_100k_R2.fastq --pairedparams 7:5:10 --outputprefix ./trial73  --h5score 0.9 --h3score 0.9 --filtering 1:50:150:20:15
