require 'spec_helper'

describe Bio::Rna do
  let(:rna) { Bio::Rna.new } 

  it 'translates rna into aminoacids' do
    sample = 'AUGGCCAUGGCGCCCAGAACUGAGAUCAAUAGUACCCGUAUUAACGGGUGA'
    result = rna.translate(sample)
    expect(result).to eq('MAMAPRTEINSTRING')
  end

  it 'finds the dna sequences that encode a given peptide' do
    sample = 'ATGGCCATGGCCCCCAGAACTGAGATCAATAGTACCCGTATTAACGGGTGA'
    peptide = 'MA'
    result = rna.encode(sample, peptide)
    expect(result).to eq(%w(ATGGCC GGCCAT ATGGCC))
  end
end
