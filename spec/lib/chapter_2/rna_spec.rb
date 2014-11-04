require 'spec_helper'

describe Bio::Rna do
  let(:rna) { Bio::Rna.new }

  it 'should translate rna into aminoacids' do
    sample = 'AUGGCCAUGGCGCCCAGAACUGAGAUCAAUAGUACCCGUAUUAACGGGUGA'
    result = rna.translate(sample)    
    expect(result).to eq('MAMAPRTEINSTRING')
  end
end
