from typing import Optional, Tuple, Union, Dict, List


class FastaParser:

    def __init__(self, file: str):
        self.__fasta = open(file, 'r')

    def __enter__(self):
        return self

    def __exit__(self, exc_type, exc_val, exc_tb):
        self.close()
        return

    def __iter__(self):
        self.__fasta.seek(0)
        return self

    def __next__(self):
        r = self.next()
        if r:
            return r
        else:  # r is None
            raise StopIteration

    def next(self) -> Optional[Tuple[str, str]]:
        """
        Returns the next read of the fasta file
        If it reaches the end of the file, return None
        """
        header = self.__fasta.readline().rstrip()[1:]
        if header == '':
            return None

        seq = ''
        while True:
            pos = self.__fasta.tell()
            line = self.__fasta.readline().rstrip()
            if line.startswith('>'):
                self.__fasta.seek(pos)
                return header, seq
            if line == '':
                return header, seq
            seq = seq + line

    def close(self):
        self.__fasta.close()


class FastaWriter:

    def __init__(self, file: str, mode: str = 'w'):
        """
        Args:
            file: path-like

            mode: 'w' for write or 'a' for append
        """
        assert mode in ['w', 'a']

        self.__fasta = open(file, mode)

    def __enter__(self):
        return self

    def __exit__(self, exc_type, exc_val, exc_tb):
        self.close()

    def __wrap(self, seq: str, length: int) -> str:
        """
        Wraps a single line of string by a specified length

        Args:
            seq: a single line of DNA or protein sequence without '\n'

            length: length of each line
        """
        if len(seq) <= length:
            return seq  # no need to wrap

        w = length
        list_ = []
        for i in range(int(len(seq)/w) + 1):
            list_.append(seq[i*w:(i+1)*w])

        return '\n'.join(list_)

    def write(self, header: str, sequence: str, wrap: int = 80):
        """
        Args:
            header: Fasta header

            sequence: DNA or protein sequence without '\n'

            wrap: length of each wrapped line for the sequence
        """
        seq = self.__wrap(sequence, wrap)
        self.__fasta.write(f'>{header}\n{seq}\n')

    def close(self):
        self.__fasta.close()
