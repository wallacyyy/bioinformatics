require 'spec_helper'

describe Bio::Rna do
  let(:peptide) { Bio::Peptide.new } 

  it 'generates the linear spectrum from a cylic peptide' do
    sample = 'NQEL'
    output = '0 113 114 128 129 242 242 257 370 371 484'
    result = peptide.linear_spectrum(sample)
    expect(result.join(' ')).to eq(output)
  end

  context 'finds the sequence of a cyclopeptide' do
    it 'runs with a simple spectrum' do
      input = peptide.format_input('0 113 128 186 241 299 314 427')
      output = %w(186-128-113 186-113-128 128-186-113 128-113-186 113-186-128 113-128-186)
      result = peptide.cyclopeptides_sequence(input)
      expect(result.sort).to eq(output.sort)
    end

    xit 'runs with a complex spectrum' do
      input = IO.read('./spec/fixtures/cyclo-input.txt').chomp 
      output = IO.read('./spec/fixtures/cyclo-output.txt').chomp 
      result = peptide.cyclopeptides_sequence(peptide.format_input(input))
      expect(result.sort).to eq(output.split.sort)
    end
  end
end
