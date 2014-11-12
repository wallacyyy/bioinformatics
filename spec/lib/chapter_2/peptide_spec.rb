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

    it 'runs with a complex spectrum' do
      input = IO.read('./spec/fixtures/cyclo-input.txt').chomp 
      output = IO.read('./spec/fixtures/cyclo-output.txt').chomp 
      result = peptide.cyclopeptides_sequence(peptide.format_input(input))
      expect(result.sort).to eq(output.split.sort)
    end
  end

  it 'calculates the score between a linear spectrum of a peptide and another' do
    sample = 'GAQLPPGWYQCGDWGIQDAKAFPHKNTMCENKTDFVQGAWDHTPM'
    input = IO.read('./spec/fixtures/score.txt').chomp 
    spectrum = peptide.format_input(input)
    expect(peptide.score(sample, spectrum)).to eq(668)
  end

  it 'finds the sequence of a cyclopeptide based on a leaderboard' do
    input = peptide.format_input('0 71 113 129 147 200 218 260 313 331 347 389 460')
    output = '129-113-147-71'
    result = peptide.leaderboard_sequence(input, 10)
    expect(result).to eq(output)
  end

  context 'trims a leaderboard' do
    it 'runs with a simple dataset' do
      leaderboard = ['LAST', 'ALST', 'TLLT', 'TQAS']
      spectrum = peptide.format_input('0 71 87 101 113 158 184 188 259 271 372')
      output = ['LAST', 'ALST']
      result = peptide.trim(leaderboard, spectrum, 2)
      expect(result).to eq(output)
    end

    it 'runs with a complex dataset' do
      leaderboard = IO.read('./spec/fixtures/trim-leaderboard.txt').chomp.split
      raw = IO.read('./spec/fixtures/trim-leaderboard-spectrum.txt').chomp
      spectrum = peptide.format_input(raw)
      output = IO.read('./spec/fixtures/trim-leaderboard-output.txt').chomp
      result = peptide.trim(leaderboard, spectrum, 5)
      expect(result.sort).to eq(output.split.sort)
    end
  end
end
