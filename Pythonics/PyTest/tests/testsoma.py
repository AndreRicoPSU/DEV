import unittest
from src.main import soma


class TesteSoma(unittest.TestCase):
    def test_retorno_soma_10_10(self):
        self.assertEqual(soma(10, 10), 20)

    def test_retorno_soma_10_20(self):
        self.assertEqual(soma(10, 10), 21)

# para rodar: python -m unittest discover -v
# video https://www.youtube.com/watch?v=ZQ60SJDACuc
