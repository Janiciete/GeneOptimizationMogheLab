# GeneOptimizationMogheLab

This project builds a plant protein mapping pipeline that converts messy gene/protein entries into structured, confidence-scored datasets linked to UniProt metadata and Phytozome protein sequences.

## Project Goal

The goal is to create a reliable plant protein mapping dataset that can later support downstream computational biology or machine learning work, such as adapting MTGNN-style approaches to plant protein research. This project does **not** train a model yet; it focuses only on building the dataset foundation.

## What the Pipeline Does

The pipeline takes a raw list of gene/protein names and:

1. Standardizes species information across databases.
2. Classifies each entry by identifier type.
3. Builds searchable UniProt and Phytozome indexes.
4. Matches entries to UniProt metadata and/or Phytozome sequences.
5. Assigns confidence and ambiguity categories.
6. Produces final confidence, coverage, and unresolved datasets.
