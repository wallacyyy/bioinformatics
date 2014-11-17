require 'spec_helper'

describe Bio::Convolution do
  let(:convolution) { Bio::Convolution.new } 
  let(:peptide)     { Bio::Peptide.new } 

  context 'calculates the spectral convolution of a spectrum' do
    it 'with a simple dataset' do
      raw = '0 137 186 323'
      spectrum = peptide.format_input(raw)
      result = convolution.spectral_convolution(spectrum)
      expect(result).to eq([137, 137, 186, 186, 323, 49].sort)
    end

    it 'with a complex dataset' do
      raw_input = IO.read('./spec/fixtures/convolute-input.txt').chomp
      output = IO.read('./spec/fixtures/convolute-output.txt').chomp
      spectrum = peptide.format_input(raw_input)
      result = convolution.spectral_convolution(spectrum.sort)
      expect(result).to eq(output.split.collect{|i| i.to_i}.sort)
    end
  end

  context 'calculates the spectral convolution of a spectrum using leaderboard' do
    it 'with a simple dataset' do
      input = peptide.format_input('57 57 71 99 129 137 170 186 194 208 228 ' +
                                   '265 285 299 307 323 356 364 394 422 493')
      m = 20
      n = 60
      output = '99-71-137-57-72-57'
      result = convolution.convolution_cyclopeptide_sequence(input, m, n)
      expect(result).to eq(output)
    end

    it 'with a complex dataset' do
      raw = IO.read('./spec/fixtures/convolute-leader.txt').chomp
      input = peptide.format_input(raw).sort
      m = 20
      n = 373
      output = '57-57-147-129-129-131-163-97-128-114-115-113-129'
      result = convolution.convolution_cyclopeptide_sequence(input, m, n)
      expect(result).to eq(output)
    end
  end
end
