#!/uufs/chpc.utah.edu/sys/installdir/julia/0.6.1/bin/julia

#= 
This script spikes a variant into a vcf file for testing
with VVP and VAAST3.

Use a bedfile to input the variant or multiple variants. 

Usage: SpikeVCF.jl myVCFfile.vcf myBedfile.bed id_index

Written: 20180510 WSW GPL
=#

if length(ARGS) < 3
    print("\n\nUsage: spikeVCF.jl myVCFfile.vcf myBedfile.bed <id_idx>

    Where the bed file contains one variant in a tab-delimited line:

    7     150649881   G    A    R    C    KCNH2   ENST00000430723

    and id_indices are the 1-based vcf sample indices to have a 
    0/1 genotype (e.g. 2,5,9).  Output is to STDOUT.\n\n")
    exit()
end

vcf = ARGS[1]
bed = ARGS[2]

s=readstring(pipeline(`grep '^#C' "$vcf"`))
ss=split(s, '\t')
num=length(ss) - 9    #number of samples

#print("Processing VCF file: $vcf\n","Found $num samples in vcf\n", "Spiking variant in $bed\n")


genos=fill("0/0", num)
mutidxs=split(ARGS[3], ',')
ac = length(mutidxs)

for i in 1:length(mutidxs)
  smutidx=shift!(mutidxs)
  mutidx=parse(Int,smutidx)
  genos[mutidx] = "0/1"
end

sgenos = join(genos, '\t')

b = []
open(ARGS[2]) do bed   
  for i in eachline(bed)
    if ismatch(r"^[X|Y|M]", i)
      println("Sorry, only autosomes are allowed in this script.")
      exit()
    end
    push!(b,i)
    println(b)
  end
end

bb = split(b[1], '\t')
cp = bb[1] * ":" *  bb[2]
alt=bb[4]
gn=bb[7]
tid=bb[8]
aa = bb[5] * "/" * bb[6]
tidc=chomp(tid)
CSQ="AC=$ac;CSQ=$alt|||$gn|||$tidc|||||||||$aa||||||"

spike = join([bb[1], bb[2], cp, bb[3], alt, "9999.99", "PASS", "$CSQ", "GT", sgenos], '\t') 

chrbed = parse(Int,bb[1]) 
posbed = parse(Int,bb[2])

open(ARGS[1]) do vcf    #read bedfile and make lookup and set

  k=0

  for i in eachline(vcf)

    if ismatch(r"^\#", i)
      println(i)
      continue
    end

    if ismatch(r"^X", i)
      println(i)
      continue
    end

    if ismatch(r"^Y", i)
      println(i)
      continue
    end


    a = split(i, '\t')
    chr = parse(Int,a[1])

    if chrbed != chr
      println(i)
      continue
    end

    if chrbed == chr
      pos = parse(Int,a[2])

      if pos < posbed 
         println(i) 
      end

      if posbed < pos
         k = k + 1; 
         if k == 1
            println(spike)
         else 
            println(i)
         end
      end
    end
  end
end
