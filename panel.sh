a=()

for x in `cat gene.txt`
do
  a+=("$x")
done

b=()

for x in `cat region.txt`
do
  b+=("$x")
done

for i in "${!b[@]}"
do
  gene="${a[$i]}"
  region="${b[$i]}"
  chrom=$(echo $region | cut -d ":" -f 1)
  fuc tabix-slice "/Users/sbslee/Desktop/bochet.gcc.biostat.washington.edu/beagle/1000_Genomes_phase3_v5a/b37.vcf/chr$chrom.1kg.phase3.v5a.vcf.gz" "$region" | fuc fuc-bgzip > "/Users/sbslee/Desktop/$gene.vcf.gz"
done
