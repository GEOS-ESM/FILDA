MUTE = 0
def printf(message, level=1, type='info', prefix = ''):
	message = f" - {prefix} {message}"
	if level != MUTE:
		print(message)
		
if __name__ == '__main__':
	print('')