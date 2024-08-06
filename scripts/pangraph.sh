#! /usr/bin/bash

#python data_utils.py


#docker pull neherlab/pangraph:latest

fasta_path=/Users/williamrshoemaker/GitHub/gutsyphage/data/vOTU-000085.fna
json_path=/Users/williamrshoemaker/GitHub/gutsyphage/data/vOTU-000085-pangraph.json


#workdir=/Users/williamrshoemaker/GitHub/gutsyphage/scripts

#docker run —rm -it —volume="$(pwd)" --user="$(id -u):$(id -g)" —workdir=${workdir} neherlab/pangraph:latest bash -c "pangraph build —circular —alpha 0 \--beta 0 ${fasta_path} > ${json_path}


#docker run --rm -it --volume="$(pwd):/workdir" \
#    --user="$(id -u):$(id -g)" --workdir=/workdir \
#    neherlab/pangraph:latest bash -c "pangraph build --circular --alpha 0 \
#    --beta 0 /workdir/vOTU-000085.fna > /workdir/vOTU-000085.json"


#This command took 5min on my laptop, which is fairly high-powered 

## RUN PANGRAPH POLISH
docker run --rm -it --volume="$(pwd):/workdir" --user="$(id -u):$(id -g)" \
    --workdir=/workdir neherlab/pangraph:latest \
    bash -c "pangraph polish /workdir/vOTU-000085-pangraph.json > /workdir/vOTU-000085-polished-pangraph.json"


#This command took 35min on my laptop