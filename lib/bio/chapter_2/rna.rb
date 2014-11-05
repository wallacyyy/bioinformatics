module Bio
  class Rna
    MASS_TABLE =
      { 'G' => 57,  'A' => 71,  'S' => 87,  'P' => 97,  'V' => 99,
        'T' => 101, 'C' => 103, 'I' => 113, 'L' => 113, 'N' => 114,
        'D' => 115, 'K' => 128, 'Q' => 128, 'E' => 129, 'M' => 131,
        'H' => 137, 'F' => 147, 'R' => 156, 'Y' => 163, 'W' => 186 }

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


    def prefix_mass_value(peptide)
      prefix_mass = [0]
      for i in 1..peptide.length do
        prefix_mass.push(prefix_mass[i - 1] + MASS_TABLE[peptide[i - 1]])
      end
      prefix_mass
    end

    def cyclic_spectrum(peptide)
      spectrum = [0]
      prefix_mass = prefix_mass_value(peptide)
      peptide_mass = prefix_mass.last

      for i in 0..peptide.length - 1 do
        for j in i + 1..peptide.length do
          spectrum.push(prefix_mass[j] - prefix_mass[i])
          if (i > 0 and j < peptide.length)
            spectrum.push(peptide_mass - (prefix_mass[j] - prefix_mass[i]))
          end
        end
      end
      spectrum.sort
    end

    def linear_spectrum(peptide)
      spectrum = [0]
      prefix_mass = prefix_mass_value(peptide)

      for i in 0..peptide.length - 1 do
        for j in i + 1..peptide.length do
          spectrum.push(prefix_mass[j] - prefix_mass[i])
        end
      end
      spectrum.sort
    end

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
