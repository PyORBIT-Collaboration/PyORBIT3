project('orbit', ['cpp'], version: '0.0.1')
inc = [
'/Applications/Xcode.app/Contents/Developer/Library/Frameworks/Python3.framework/Versions/3.9/Headers',
'/opt/homebrew/include'
]

link_a = ['-L/Applications/Xcode.app/Contents/Developer/Library/Frameworks/Python3.framework/Versions/3.9/lib',
      '-L/opt/homebrew/lib',
      '-lpython3.9', '-lfftw3'
]

subdir('src')
subdir('py/orbit')