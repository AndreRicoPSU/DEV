from word_extractor import BaseExtractor


if __name__ == "__main__":
	
	# url for word list (huge)
	url = 'https://raw.githubusercontent.com/dwyl/english-words/master/words.txt'
	
	# word list in array
	array = ['one', 'two', 'three', 'four', 'five']
	
	# initializing BaseExtractor object
	extractor = BaseExtractor()
	
	# calling parse_url function
	extractor.parse_url(url)
	
	# calling pasrse_list function
	extractor.parse_list(array)
