#! //usr/bin/ruby
arg = ARGV.join("\s")
arg ="#{arg}\s"

/-i\s/ =~ arg
prefix = $'
/\s/ =~ prefix
sam = $`

`rm -r sam_demutiplexed`
`mkdir sam_demutiplexed`
fh = File.open("#{sam}","r")
while !fh.eof
	line = fh.gets.chomp
	barcode = line.split("\t")[0].split("_")[1]
	`echo "#{line}" >> sam_demutiplexed/#{barcode}.sam`
end