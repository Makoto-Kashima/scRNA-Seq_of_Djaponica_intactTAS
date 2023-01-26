#!/usr/bin/ruby
arg = ARGV.join("\s")
arg ="#{arg}\s"
puts arg
#set  the list of seqs
/-i\s/ =~ arg
inputs = $'
/\s/ =~ inputs
inputs = $`

#set  refference
/-d\s/ =~ arg
db = $'
/\s/ =~ db
db = $`

#set  refference
/-a\s/ =~ arg
cpu = $'
/\s/ =~ cpu
cpu = $`

`mkdir bwa_mem`
out = File.open("bwa_mem.log","a")
fh = File.open("#{inputs}","r")
while !fh.eof
input = fh.gets.chomp
puts input
	out.puts "bwa mem -t #{cpu} #{db} #{input} | samtools view -@ #{cpu} -b  | samtools sort -@ #{cpu} -m  2G > bwa_mem/#{input.split("/")[1].gsub("fq.gz","bam")}"
	`bwa mem -t #{cpu} #{db} #{input} | samtools view -@ #{cpu} -b  | samtools sort -@ #{cpu} -m  2G > bwa_mem/#{input.split("/")[1].gsub("fq.gz","bam")}`
	`samtools index bwa_mem/#{input.split("/")[1].gsub("fq.gz","bam")}`
end
out.close
fh.close
