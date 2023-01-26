#! /usr/bin/ruby
arg = ARGV.join("\s")
arg ="#{arg}\s"
#set  the list of seqs
/-i\s/ =~ arg
input = $'
/\s/ =~ input
list = $`

#set database name
/-f\s/ =~ arg
refseq = $'
/\s/ =~ refseq	
refseq = $`

fh = File.open("#{list}","r")
`mkdir bwa_mem_split`
while !fh.eof
	input = fh.gets.chomp
	input2 = input.split("/")[1]
	/.sam/ =~ input2
	input2 = $`
	samtools = `samtools view -@ 60 -Shb #{input} -T #{refseq} -o bwa_mem_split/#{input2}.bam`
	samtools = `samtools index bwa_mem_split/#{input2}.bam`
end