#!/bin/bash

## Test dataset - polyA RNA
../bin/RNApip_cluster -i SAMPLE_PLAN_B249 -o /data/tmp/RNAPIP_OP_B249/ -c CONFIG_Human -n testop_B249

## Test dataset - total RNA
../bin/RNApip_cluster -i SAMPLE_PLAN_L81 -o /data/tmp/RNAPIP_OP_L81/ -c CONFIG_Human -n testop_L81

## Test dataset - Mouse samples
../bin/RNApip_cluster -i SAMPLE_PLAN_Gendrel -o /data/tmp/RNAPIP_Gendrel -c CONFIG_Mouse -n testop_Gendrel

