require 'spec_helper'

xdescribe Bio::DnaComplement do
  let(:reverse_input)       { IO.read('./spec/fixtures/reverse-input.txt') }
  let(:reverse_output)      { IO.read('./spec/fixtures/reverse-output.txt') }
  let(:clump_output)        { IO.read('./spec/fixtures/reverse-output.txt') }

  describe 'gives the reverse complement' do
    it 'runs with a simple dna string' do   
      complementer = Bio::DnaComplement.new('AAAACCCGGT')
      expect(complementer.reverse).to eq('ACCGGGTTTT')
    end

    it 'runs with a greater dna string' do
      complementer = Bio::DnaComplement.new(reverse_input)
      expect(complementer.reverse).to eq(reverse_output.chomp)
    end
  end

  describe 'find the starting position of a pattern' do
    it 'runs with a simple dna string' do
      complementer = Bio::DnaComplement.new('GATATATGCATATACTT')
      expect(complementer.positions('ATAT')).to eq([1, 3, 9])
    end
  end
end
