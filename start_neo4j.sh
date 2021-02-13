#!/bin/bash

docker run \
    -p 0.0.0.0:8474:7474 \
    -p 0.0.0.0:8473:7473 \
    -p 0.0.0.0:8687:7687 \
    --name ensep_neo4j \
    -v $PWD/data:/data \
    --env NEO4J_AUTH=neo4j/ensep \
    --env NEO4J_dbms_memory_pagecache_size=4G \
    --user="$(id -u):$(id -g)" \
    neo4j
