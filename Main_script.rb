require 'bio'
require 'rest-client'

puts ("Creating the gff3 file...")
out = File.open("Allrepeats.gff3", "w")
out.puts ("##gff-version 3")
out.close

puts ("Creating the file with genes with no repeats...")
out = File.open("No_repeats_CTTCTT.txt", "w")
out.puts("This file contains a list with the genes that doesn't have any CTTCTT repetition:")
out.close

puts ("Creating the gff3 file with the chr coordinates...")
out = File.open("Allrepeatschr.gff3", "w")
out.puts("##gff-version 3")
out.close

def fetch(url, headers = {accept: "*/*"}, user = "", pass="")
  response = RestClient::Request.execute({
    method: :get,
    url: url.to_s,
    user: user,
    password: pass,
    headers: headers})
  return response
  
  rescue RestClient::ExceptionWithResponse => e
    $stderr.puts e.inspect
    response = false
    return response  # now we are returning 'False', and we will check that with an \"if\" statement in our main code
  rescue RestClient::Exception => e
    $stderr.puts e.inspect
    response = false
    return response  # now we are returning 'False', and we will check that with an \"if\" statement in our main code
  rescue Exception => e
    $stderr.puts e.inspect
    response = false
    return response  # now we are returning 'False', and we will check that with an \"if\" statement in our main code
end

def create_feature(position,strand)
  if strand == "+"
  newfeat = Bio::Feature.new("repeat", "#{position}")
  newfeat.append(Bio::Feature::Qualifier.new('repeat motif','CTTCTT'))
  newfeat.append(Bio::Feature::Qualifier.new('strand', "#{strand}"))      
  elsif strand == "-"
  newfeat = Bio::Feature.new("repeat", "#{position}")
  newfeat.append(Bio::Feature::Qualifier.new('repeat motif','CTTCTT'))
  newfeat.append(Bio::Feature::Qualifier.new('strand', "#{strand}"))      
  end 
  return newfeat 
end

def get_genes(gene_id)
  
  gene_embl = {}
  
  response = fetch("http://www.ebi.ac.uk/Tools/dbfetch/dbfetch?db=ensemblgenomesgene&format=embl&id=#{gene_id}")
          
    record = Bio::EMBL.new(response.body)
    record = record.to_biosequence
    gene_embl[gene_id] = record
         
  return gene_embl
end

def pattern(pattern)
  seq = Bio::Sequence.auto(pattern)
  search = Regexp.new(seq.to_re)
  return search
end

def write_gffile(file_name,cromnum,featurename,position,strand,attr)
  
  out = File.open(file_name, "a")
  out.puts("Chr#{cromnum}\tEMBL\t#{featurename}\t#{(position.to_i) + 1}\t#{(position.to_i) + 6}\t.\t#{strand}\t.\t#{attr}")
  out.close  
end

def write_gffilechr(file_name,cromnum,featurename,position,strand,attr,chr)
  
  out = File.open(file_name, "a")
  out.puts("#{cromnum}\tEMBL\t#{featurename}\t#{(position.to_i) + 1 + chr.to_i}\t#{(position.to_i) + 6 + chr.to_i}\t.\t#{strand}\t.\t#{attr}")
  out.close  
end

def create_outfile(file_name,gene_id)
  
  out = File.open(file_name, "a")
  out.puts("#{gene_id}")
  out.close  
end

def analyse_exons(hash)
  new_feats = []
  match_pos = []
  match_complpos = []  
hash.each_value do |gene|
  exon_positions = []
  exon_complpositions = []  
  cromnum = gene.entry_id
  chromosome = gene.accessions[0].split(":")
  chr = chromosome[3]
  gene.features.each do |feature|
    if feature.feature == "exon"
      if /:/.match(feature.position)
        next
      elsif /complement/.match(feature.position)
        exon_complpositions.append(feature.position)
      else exon_positions.append(feature.position)
      end 
    end    
  end
  exon_positions.each do |position|
    position = position.split("..")
    position[0] = (position[0].to_i) -1
    positiontosum = position[0]
    position[1] = (position[1].to_i) -1
    position = position[0]..position[1]    
    match = gene.seq[position].gsub(pattern("CTTCTT")).map{Regexp.last_match.begin(0)}
    match.each do |match|
      match = match + positiontosum
      match_pos.append(match)
    end 
   end
  match_pos = match_pos.uniq   
  match_pos.each do |position| 
    newfeat = create_feature(position,"+")
    new_feats.append(newfeat)
  end 
  exon_complpositions.each do |position|
    position = position.split("(")
    position.delete(position[0])
    position = position[0].split(")")
    position = position[0].split("..")
    position[0] = (position[0].to_i) -1
    positiontosum = position[0]
    position[1] = (position[1].to_i) -1
    position = position[0]..position[1]
    match = gene.seq[position].gsub(pattern("AAGAAG")).map{Regexp.last_match.begin(0)}
    match.each do |match|
      match = match + positiontosum
      match_complpos.append(match)
    end
  end
  match_complpos = match_complpos.uniq
  match_complpos.each do |position| #unless match_pos == ""
    newfeat = create_feature(position,"-")
    new_feats.append(newfeat)
  end 
  new_feats.each do |feature|
    if feature.feature == "repeat"
      qual = feature.assoc
      position = feature.position
      featurename = feature.feature
      strand = qual['strand']
      attr = qual['repeat motif']
      write_gffile("Allrepeats.gff3",cromnum,featurename,position,strand,attr)
      write_gffilechr("Allrepeatschr.gff3",cromnum,featurename,position,strand,attr,chr)
    end
  end   
end
return match_pos.length + match_complpos.length
end 

File.open("ArabidopsisSubNetwork_GeneList.txt", "r").each do |gene_id|
  gene_id = gene_id.strip
  genes_embl = get_genes(gene_id)
  matchnum = analyse_exons(genes_embl)
  
  if matchnum == 0
    create_outfile("No_repeats_CTTCTT.txt",gene_id)    #code
  end
  
end
puts("Done")
