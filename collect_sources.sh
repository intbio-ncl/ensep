#!/bin/bash
# Download the data files.
wget https://nextcloud.bioswarm.net/s/cggB45agtonZPgA/download
# Unzip the datafiles.
unzip download
# Remove the zip
rm download

# Change directory
cd datafiles
# Decompress stuff
# MetaCyc
tar -zxf metacyc.tar.gz
rm metacyc.tar.gz
# KEGG
tar -zxf kegg.tar.gz
rm kegg.tar.gz
# ChEBI
tar -zxf chebi.tar.gz
rm chebi.tar.gz
# Rhea
tar -zxf rhea.tar.gz
rm rhea.tar.gz