require 'spec_helper'

describe Bio::Rna do
  it 'translates rna into aminoacids' do
    sample = 'AUGGCCAUGGCGCCCAGAACUGAGAUCAAUAGUACCCGUAUUAACGGGUGA'
    result = Bio::Rna.new.translate(sample)
    expect(result).to eq('MAMAPRTEINSTRING')
  end

  it 'finds the dna sequences that encode a given peptide' do
    sample = 'ATGGCCATGGCCCCCAGAACTGAGATCAATAGTACCCGTATTAACGGGTGA'
    peptide = 'MA'
    result = Bio::Rna.new.encode(sample, peptide)
    expect(result).to eq(%w(ATGGCC GGCCAT ATGGCC))
  end

  it 'generates the linear spectrum from a cylic peptide' do
    peptide = 'NQEL'
    output = '0 113 114 128 129 242 242 257 370 371 484'
    result = Bio::Rna.new.linear_spectrum(peptide)
    expect(result.join(' ')).to eq(output)
  end

end
