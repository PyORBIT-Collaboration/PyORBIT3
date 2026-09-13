#ifndef __STRING_UTILS_H_
#define __STRING_UTILS_H_

#include <string>
#include <vector>


//=======================================
//Usage
// int nTokens = StringUtils::Tokenize("....",vector<string>& tokens, ", :<>");
//or for space-delimiter
//int nTokens = StringUtils::Tokenize("....",vector<string>& tokens);
//=======================================

/** The string utilits collection. */

namespace OrbitUtils{
	namespace StringUtils
	{
		/**
		  The Tokenizer. It returns the vector with tokens.
		  The delimiter is optional, and by default it is a space.
		*/
		int Tokenize(const std::string& str,std::vector<std::string>& tokens,  const std::string& delimiters = " ");
	};
};
#endif
