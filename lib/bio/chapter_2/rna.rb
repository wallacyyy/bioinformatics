module Bio
  class Rna
    GENETIC_CODE = 
      { 'AAA' => 'K', 'AAC' => 'N', 'AAG' => 'K', 'AAU' => 'N', 
        'ACA' => 'T', 'ACC' => 'T', 'ACG' => 'T', 'ACU' => 'T', 
        'AGA' => 'R', 'AGC' => 'S', 'AGG' => 'R', 'AGU' => 'S', 
        'AUA' => 'I', 'AUC' => 'I', 'AUG' => 'M', 'AUU' => 'I', 
        'CAA' => 'Q', 'CAC' => 'H', 'CAG' => 'Q', 'CAU' => 'H', 
        'CCA' => 'P', 'CCC' => 'P', 'CCG' => 'P', 'CCU' => 'P', 
        'CGA' => 'R', 'CGC' => 'R', 'CGG' => 'R', 'CGU' => 'R', 
        'CUA' => 'L', 'CUC' => 'L', 'CUG' => 'L', 'CUU' => 'L', 
        'GAA' => 'E', 'GAC' => 'D', 'GAG' => 'E', 'GAU' => 'D', 
        'GCA' => 'A', 'GCC' => 'A', 'GCG' => 'A', 'GCU' => 'A', 
        'GGA' => 'G', 'GGC' => 'G', 'GGG' => 'G', 'GGU' => 'G', 
        'GUA' => 'V', 'GUC' => 'V', 'GUG' => 'V', 'GUU' => 'V', 
        'UAA' => '',  'UAC' => 'Y', 'UAG' => '',  'UAU' => 'Y', 
        'UCA' => 'S', 'UCC' => 'S', 'UCG' => 'S', 'UCU' => 'S', 
        'UGA' => '',  'UGC' => 'C', 'UGG' => 'W', 'UGU' => 'C', 
        'UUA' => 'L', 'UUC' => 'F', 'UUG' => 'L', 'UUU' => 'F' }

    def translate(sample)
      aminoacids = ''
      i = 0

      while (i + 2 < sample.length) do
        aminoacids << GENETIC_CODE[sample[i..i + 2]]
        i += 3 
      end
      aminoacids
    end

    def encode(sample, peptide)
      dna = DnaComplement.new
      patterns = []
      size = peptide.length
      i = 0

      while (i + size * 3 - 1 < sample.length) do
        sequence = sample[i..i + size * 3 - 1]
        reversed = dna.reverse(sequence)
        if (translate(dna.to_rna(sequence)) == peptide) or 
           (translate(dna.to_rna(reversed)) == peptide)
          patterns.push(sequence)
        end
        i += 1
      end
      patterns
    end
  end
end
