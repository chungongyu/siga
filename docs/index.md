# Home

[![Build Status](https://travis-ci.org/chungongyu/siga.svg?branch=master)](https://travis-ci.org/chungongyu/siga)
![GitHub release](https://img.shields.io/github/release/chungongyu/siga.svg)
[![GitHub license](https://img.shields.io/github/license/chungongyu/siga.svg)](https://github.com/chungongyu/siga)
[![GitHub issues](https://img.shields.io/github/issues/chungongyu/siga.svg)](https://github.com/chungongyu/siga/issues)

## Overview
SIGA is an open-source _de novo_ assembly toolkit containing various assembly pipelines which are totally compatible with [SGA](https://github.com/jts/sga) file format. Similar to [SGA](https://github.com/jts/sga), it is designed as a [modular set of programs](commands.html), which are used to form an assembly pipeline. A description of the SIGA design is found [here](design.html). 

## Data exploration and quality control

It is highly recommended that you run the 'preqc' module on your data prior to assembly. This module will give you information about the read quality and genome characteristics. See [this page](https://github.com/chungongyu/siga/wiki/preqc) for more information.

## First steps

The source directory contains [examples](https://github.com/chungongyu/siga/tree/master/examples) of real assemblies using SIGA. You should read these scripts or (better) download the data for one of the smaller genomes (I recommend the C. elegans data set) and run the example yourself. This will help you get understand the SIGA pipeline so you can run the assembler effectively on your own data.

## FAQs

It is highly recommended that you read the [SIGA FAQ page](https://github.com/chungongyu/siga/wiki/FAQ).

## Further help

There is a mailing list for siga on [google groups](http://groups.google.com/group/siga-users). 
