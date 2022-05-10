import re
from collections import defaultdict
from itertools import combinations

from mrjob.job import MRJob
from mrjob.step import MRStep

WORD_RE = re.compile(r"[\w']+")


class MRRelativeFreq(MRJob):
    denoms = defaultdict(int)

    def steps(self):
        return [
            MRStep(mapper=self.mapper, combiner=self.combiner, reducer=self.reducer),
            # MRStep(reducer=self.reducer_s2),
        ]

    def mapper(self, _, line):
        words = WORD_RE.findall(line)
        for (x, y) in combinations(words, 2):
            if x != y:
                yield ((x.lower(), "*"), 1)
                yield ((x.lower(), y.lower()), 1)

    def combiner(self, pair, counts):
        yield (pair, sum(counts))

    def reducer(self, pair, counts):
        count = sum(counts)
        x, y = pair
        if y == "*":
            self.denoms[x] = count
        else:
            yield ((x, y), count)

    def reducer_s2(self, pair, ycnt):
        x, y = pair
        lkup = self.denoms[x]
        yield (pair, round((sum(ycnt) / lkup), 2))


if __name__ == "__main__":
    MRRelativeFreq.run()
