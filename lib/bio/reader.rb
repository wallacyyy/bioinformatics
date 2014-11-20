module Bio
  class Reader
    def initialize(path)
      @path = path
    end

    def read_lines
      lines = []
      IO.foreach(@path) do |line|
        lines.push(line.chomp)
      end
      lines
    end
  end
end
