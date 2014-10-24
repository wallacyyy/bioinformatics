require 'spec_helper'

describe Bio::DnaComplement do
  let(:reverse_input_path)  { './spec/fixtures/reverse-input.txt' }
  let(:reverse_input)       { IO.readlines(reverse_input_path).first }
  let(:reverse_output_path) { './spec/fixtures/reverse-output.txt' }
  let(:reverse_output)      { IO.readlines(reverse_output_path).first }
  let(:clump_path)          { './spec/fixtures/reverse-output.txt' }
  let(:clump_output)        { IO.readlines(reverse_output_path).first }


  describe 'gives the reverse complement' do
    it 'runs with a simple dna string' do   
      complementer = Bio::DnaComplement.new('AAAACCCGGT')
      expect(complementer.reverse).to eq('ACCGGGTTTT')
    end

    it 'runs with a greater dna string' do
      complementer = Bio::DnaComplement.new(reverse_input)
      expect(complementer.reverse).to eq(reverse_output.gsub("\n", ""))
    end
  end

  describe 'find the starting position of a pattern' do
    it 'runs with a simple dna string' do
      complementer = Bio::DnaComplement.new('GATATATGCATATACTT')
      expect(complementer.positions('ATAT')).to eq([1, 3, 9])
    end
  end

  it 'find oriC clumps' do
    complementer = Bio::DnaComplement.new('CGGACTCGACAGATGTGAAGAACGACAATGTG' +
                                          'AAGACTCGACACGACAGAGTGAAGAGAAGAGGAAACATTGTAA')
    expect(complementer.clumps(5, 50, 4)).to eq(['CGACA', 'GAAGA'])
  end
end

