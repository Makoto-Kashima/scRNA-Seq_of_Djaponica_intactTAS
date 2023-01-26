#!/usr/bin/ruby
arg = ARGV.join("\s")
arg ="#{arg}\s"

/-i\s/ =~ arg
input1 = $'
/\s/ =~ input1
input1 = $`

/-f\s/ =~ arg
fasta = $'
/\s/ =~ fasta
fasta = $`

#set num of CPU 
/-a\s/ =~ arg
a = $'
/\s/ =~ a
a = $`.to_i
if a == nil
	a=64
end

`mkdir salmon`
out = File.open("salmon.log","a")
fh = File.open("#{input1}","r")


rm = `rm #{input1}_list`
thread = []
files = []
t = 4
split = `split -n l/#{t} #{input1} #{input1}_`
ls = `ls #{input1}_[a-z][a-z] > #{input1}_list`
fh = File.open("#{input1}_list", "r")
i = 0
while !fh.eof
	file = fh.gets.chomp	
	thread.push "#{file}"
	files.push "#{file}"
end
def Flexible_UMI_Demultiplexing(input, fasta, out)
	fh = File.open("#{input}","r")
	while !fh.eof
		input0 = fh.gets.chomp
		input = input0.split("/")[1]
		out.puts "salmon quant -t #{fasta} -l IU -p 16 -a #{input0} -o salmon/#{input.split(".")[0]}"
		`salmon quant -t #{fasta} -l IU -p 16 -a #{input0} -o salmon/#{input.split(".")[0]}`
		`rm -r salmon/#{input.split(".")[0]}/logs`
	end
end
t.times do |j|
	thread[j] = Thread.new do
 	file = files[j]
	puts file
	Flexible_UMI_Demultiplexing(file, fasta, out)
	end
	sleep(1)
end


t.times do |z|
 	t.times do |z|
		thread[z].join

	end
end
rm = `rm #{input1}_[a-z][a-z] #{input1}_list`
out.close
fh.close
`rm salmon/*/logs/*`

