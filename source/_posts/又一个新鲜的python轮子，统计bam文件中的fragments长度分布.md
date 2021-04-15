---
title: 又一个新鲜的python轮子，统计bam文件中的fragments长度分布
date: 2021-04-16 02:04:34
tags:
  - 原创
  - 生信
categories:
  - [python]
  - [生信入门]
description: 输入bam文件，输出两列，一列是fragment的长度，一列是数量。不愧是python，就是慢。
---

[RSeQC](https://github.com/MonashBioinformaticsPlatform/RSeQC)的[RNA_fragment_size.py](https://github.com/MonashBioinformaticsPlatform/RSeQC/blob/master/rseqc/modules/RNA_fragment_size.py)虽然可以统计fragments的长度，但只能输出每个比对到每个基因的片段的均值，中位数等信息：

> calculate fragment size for each gene/transcript. For each transcript/gene, it Will report:
> 1) # of fragment that was used.
> 2) mean of fragment size
> 3) median of fragment size
> 4) stdev of fragment size

所以我在它的基础上整了一个烂活，可以统计fragment的长度，方便生成长度分布直方图。

### 使用方法：

```shell
usage: fragment_sizes.py [-h] -i INPUT_FILE [-o OUTPUT_FILE] -r REFGENE_BED
                         [-q MAP_QUAL] [-n FRAGMENT_NUM] [-t TEMP_FILE]

calculate fragment size for each gene/transcript.

optional arguments:
  -h, --help            show this help message and exit
  -i INPUT_FILE, --input INPUT_FILE
                        Input BAM file (default: None)
  -o OUTPUT_FILE, --output OUTPUT_FILE
                        Output file (default: None)
  -r REFGENE_BED, --refgene REFGENE_BED
                        Reference gene model in BED format. Must be strandard
                        12-column BED file. [required] (default: None)
  -q MAP_QUAL, --mapq MAP_QUAL
                        Minimum mapping quality (phred scaled) for an
                        alignment to be called "uniquely mapped". default=30
                        (default: 30)
  -n FRAGMENT_NUM, --frag-num FRAGMENT_NUM
                        Minimum number of fragment. default=1 (default: 1)
  -t TEMP_FILE, --temp TEMP_FILE
                        Output file to store processing data, support *.gz.
                        default=devnull (default: /dev/null)
```

讲人话就是：

```shell
fragment_sizes.py -i $bam文件 -r $bed12 -o $输出文件 -t $缓存文件
```

输出文件-o可以不加，然后输出就会打印到标准输出。缓存文件-t实际上就是每一行一个fragment的长度，也可以不加，支持gzip压缩格式，后缀加“.gz"就行了。

### 源码：

```python
#!/usr/bin/env python3
'''
calculate fragment size for each gene/transcript.
'''

import os
import sys
import gzip
import pysam
import argparse
import pandas as pd
from collections import Counter


def overlap_length(lst1, lst2):
	l = 0
	for x in lst1:
		for y in lst2:
			l += len(range(max(x[0], y[0]), min(x[-1], y[-1]) + 1))
	return l


def fragment_size(bedfile, samfile, qcut=30, ncut=1, temp=os.devnull):
	'''calculate the fragment size for each gene'''
	for line in open(bedfile, 'r'):
		exon_range = []
		if line.startswith(('#','track','browser')):
			continue
		fields = line.split()
		chrom = fields[0]
		tx_start = int(fields[1])
		tx_end = int(fields[2])
		geneName = fields[3]
		trand = fields[5].replace(" ", "_")
		exon_starts = map(int, fields[11].rstrip(',\n').split(','))
		exon_starts = map((lambda x: x + tx_start), exon_starts)
		exon_ends = map(int, fields[10].rstrip(',\n').split( ',' ))
		exon_ends = map((lambda x, y: x + y), exon_starts, exon_ends)

		for st, end in zip(exon_starts, exon_ends):
			exon_range.append([st + 1,end + 1])
		try:
			alignedReads = samfile.fetch(chrom, tx_start, tx_end)
		except:
			continue

		frag_sizes = []
		for aligned_read in alignedReads:
			if not aligned_read.is_paired:				#skip single sequencing
				continue
			if aligned_read.is_read2:
				continue
			if aligned_read.mate_is_unmapped:
				continue
			if aligned_read.is_qcfail:
				continue                                #skip low quanlity
			if aligned_read.is_duplicate:
				continue                                #skip duplicate read
			if aligned_read.is_secondary:
				continue                                #skip non primary hit
			if aligned_read.mapq < qcut:
				continue

			read_st = aligned_read.pos
			mate_st = aligned_read.pnext
			if read_st > mate_st:
				(read_st, mate_st) = (mate_st, read_st)
			if read_st < tx_start or mate_st > tx_end:
				continue
			read_len = aligned_read.qlen
			map_range = [[read_st+1, mate_st]]
			frag_len = overlap_length(exon_range, map_range) + read_len
			frag_sizes.append(frag_len)
		if len(frag_sizes) < ncut:
			continue
		else:
			for i in frag_sizes:
				temp.write(f"{i}\n".encode())
				yield i


if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description=__doc__,
        formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument("-i", "--input", dest="input_file", type=str, required=True, help="Input BAM file")
    parser.add_argument("-o", "--output", dest="output_file", type=str, help="Output file")
    parser.add_argument("-r","--refgene", dest="refgene_bed", type=str, required=True,
        help="Reference gene model in BED format. Must be strandard 12-column BED file. [required]")
    parser.add_argument("-q","--mapq", dest="map_qual", type=int, default=30,
        help="Minimum mapping quality (phred scaled) for an alignment to be called \"uniquely mapped\". default=30")
    parser.add_argument("-n","--frag-num", dest="fragment_num", type=int, default=1,
        help="Minimum number of fragment. default=1")
    parser.add_argument("-t","--temp", dest="temp_file", type=str, default=os.devnull,
		help="Output file to store processing data, support *.gz. default=devnull")

    args = parser.parse_args()

    if not os.path.exists(args.input_file + '.bai'):
        print >>sys.stderr, "cannot find index file of input BAM file"
        print >>sys.stderr, args.input_file + '.bai' + " does not exists"
        sys.exit(0)

    temp_writer = gzip.open(args.temp_file, "wb") if args.temp_file.endswith(".gz") else open(args.temp_file, "wb")
    fragment_sizes = Counter(i for i in fragment_size(
		args.refgene_bed,
		pysam.Samfile(args.input_file),
		args.map_qual,
		args.fragment_num,
		temp_writer
		))

    frag_sizes_df = pd.DataFrame({"Length": fragment_sizes.keys(), "Count": fragment_sizes.values()})
    frag_sizes_df.sort_values(by="Length", ascending=True)

    if args.output_file:
        frag_sizes_df.to_csv(args.output_file, sep="\t", header=True, index=False)
    else:
        print("Length\tCount")
        for _, i in frag_sizes_df.iterrows():
            print(i["Length"], i["Count"], sep="\t")

```

