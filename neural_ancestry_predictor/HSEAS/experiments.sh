#!/bin/bash
source ../scripts/start_genomics_universal.sh
export CUDA_VISIBLE_DEVICES=0

python3 neural_ancestry_predictor.py --config HSEAS/k14.yaml > HSEAS/results_k14.txt 2>&1
python3 neural_ancestry_predictor.py --config HSEAS/k15.yaml > HSEAS/results_k15.txt 2>&1
# python3 neural_ancestry_predictor.py --config HSEAS/k11.yaml > HSEAS/results_k11.txt 2>&1
# python3 neural_ancestry_predictor.py --config HSEAS/k13.yaml > HSEAS/results_k13.txt 2>&1
# python3 neural_ancestry_predictor.py --config HSEAS/k21.yaml > HSEAS/results_k21.txt 2>&1
# python3 neural_ancestry_predictor.py --config HSEAS/k23.yaml > HSEAS/results_k23.txt 2>&1
# python3 neural_ancestry_predictor.py --config HSEAS/fc128.yaml > HSEAS/results_fc128.txt 2>&1
# python3 neural_ancestry_predictor.py --config HSEAS/fc512.yaml > HSEAS/results_fc512.txt 2>&1
# python3 neural_ancestry_predictor.py --config HSEAS/w8.yaml > HSEAS/results_w8.txt 2>&1
# python3 neural_ancestry_predictor.py --config HSEAS/w128.yaml > HSEAS/results_w128.txt 2>&1
# python3 neural_ancestry_predictor.py --config HSEAS/w512.yaml > HSEAS/results_w512.txt 2>&1
# python3 neural_ancestry_predictor.py --config HSEAS/9g.yaml > HSEAS/results_9g.txt 2>&1
# python3 neural_ancestry_predictor.py --config HSEAS/7g.yaml > HSEAS/results_7g.txt 2>&1
