/**
 * @file RNA_secondary.cpp
 * @brief This file contains a C++ program that computes the optimal secondary structure of an RNA sequence using dynamic programming.
 * @author Y.S Sandeep
 * @author Suman Sekhar Sahoo
 * @author Raghuram Venkatesan
 * @author Bhaarat K
 * @author Kamalesh Ram
 * @date 2024-04-28
 */

#include <bits/stdc++.h>
using namespace std;

/**
 * @brief Checks if two characters representing nucleotides are complementary (can form a base pair).
 *
 * This function takes two characters representing nucleotides (A, U, G, or C) and checks if they are complementary
 * and can form a base pair. It does this by checking if the two characters are A and U, or G and C, respectively.
 *
 * @param a The first nucleotide character.
 * @param b The second nucleotide character.
 * @return True if the two nucleotides can form a base pair, false otherwise.
 */
bool isMatching(char a, char b)
{
  unordered_set<char> st;
  st.insert(a);
  st.insert(b);
  if ((st.find('A') != st.end() && st.find('U') != st.end()) || (st.find('G') != st.end() && st.find('C') != st.end()))
  {
    return true;
  }
  return false;
}

/**
 * @brief Prints a 2D vector of integers.
 *
 * This function takes a 2D vector of integers and prints its contents to the console. It iterates over each element
 * of the vector and prints it on a new line after each row.
 *
 * @param v The 2D vector to be printed.
 */
void printDP(vector<vector<int>> v)
{
  for (int i = 0; i < v.size(); i++)
  {
    for (int j = 0; j < v.size(); j++)
    {
      cout << v[i][j] << " ";
    }
    cout << endl;
  }
}

/**
 * @brief Computes the optimal secondary structure of an RNA sequence using dynamic programming.
 *
 * This function implements the Nussinov algorithm to compute the optimal secondary structure of an RNA sequence
 * using dynamic programming. The algorithm fills a 2D table `opt` where `opt[i][j]` represents the maximum number
 * of base pairs that can be formed by the substring from index `i` to index `j` in the RNA sequence.
 *
 * The function first initializes the table for small substrings of length up to 4, where no base pairs can be formed.
 * Then, it iterates over larger substrings and computes the maximum number of base pairs that can be formed by
 * considering all possible pairs of indices `k` and `j` such that the nucleotides at those indices can form a base
 * pair. It updates the `opt` table accordingly.
 *
 * @param str The RNA sequence as a string.
 * @return The dynamic programming table representing the optimal secondary structure.
 */
vector<vector<int>> optimalMatching(string str)
{
  int n = str.size();
  vector<vector<int>> opt(n, vector<int>(n, -1));
  // Initialize the table for small substrings
  for (int diff = 0; diff <= min(4, n - 1); diff++)
  {
    for (int i = 0; i < n - diff; i++)
    {
      int j = i + diff;
      opt[i][j] = 0;
      opt[j][i] = 0;
    }
  }
  // Compute the table for larger substrings
  for (int diff = 5; diff < n; diff++)
  {
    for (int i = 0; i < n - diff; i++)
    {
      int j = i + diff;
      vector<int> matched_w_j;
      for (int t = i; t <= j - 5; t++)
      {
        if (isMatching(str[j], str[t]))
        {
          matched_w_j.push_back(t);
        }
      }
      if (matched_w_j.size())
      {
        int max_val = -1;
        for (auto k : matched_w_j)
        {
          if (k == 0)
          {
            max_val = max(max_val, 1 + opt[k + 1][j - 1]);
          }
          else
          {
            max_val = max(max_val, 1 + opt[i][k - 1] + opt[k + 1][j - 1]);
          }
        }
        opt[i][j] = max(opt[i][j - 1], max_val);
      }
      else
      {
        opt[i][j] = opt[i][j - 1];
      }
    }
  }
  return opt;
}

/**
 * @brief Performs traceback to reconstruct the optimal secondary structure from the dynamic programming table.
 *
 * This function performs a traceback on the dynamic programming table `dp` to reconstruct the optimal secondary
 * structure of the RNA sequence. It takes the starting and ending indices `i` and `j`, the RNA sequence `str`,
 * the dynamic programming table `dp`, and a reference to a vector `ans` that will store the bonded pairs of indices.
 *
 * The function works recursively by checking if the maximum number of base pairs `dp[i][j]` is the same as
 * `dp[i][j-1]`, in which case it recursively calls itself with `j-1` as the new ending index. Otherwise, it finds
 * the index `k` such that the nucleotides at indices `k` and `j` can form a base pair, and the value `dp[i][j]` is
 * the sum of 1 and the maximum number of base pairs for the substrings from `i` to `k-1` and from `k+1` to `j-1`.
 * It then adds the pair `(k, j)` to the `ans` vector and recursively calls itself with the updated indices.
 *
 * @param i The starting index of the substring.
 * @param j The ending index of the substring.
 * @param str The RNA sequence as a string.
 * @param dp The dynamic programming table.
 * @param ans A reference to a vector that will store the bonded pairs of indices.
 */
void traceback(int i, int j, string str, vector<vector<int>> dp, vector<pair<int, int>> &ans)
{
  if (j <= i)
  {
    return;
  }
  else if (dp[i][j] == dp[i][j - 1])
  {
    traceback(i, j - 1, str, dp, ans);
  }
  else
  {
    vector<int> kj_matched;
    for (int k = i; k < j; k++)
    {
      if (isMatching(str[k], str[j]))
      {
        kj_matched.push_back(k);
      }
    }
    for (auto k : kj_matched)
    {
      if (k == 0)
      {
        if (dp[i][j] == (1 + dp[k + 1][j - 1]))
        {
          ans.push_back(make_pair(0, j));
          traceback(k + 1, j - 1, str, dp, ans);
          return;
        }
      }
      else
      {
        if (dp[i][j] == (1 + dp[i][k - 1] + dp[k + 1][j - 1]))
        {
          ans.push_back(make_pair(k, j));
          traceback(i, k - 1, str, dp, ans);
          traceback(k + 1, j - 1, str, dp, ans);
          return;
        }
      }
    }
  }
}


/**
 * @brief Prints an array of pairs of integers.
 *
 * This function takes an array (vector) of pairs of integers and prints them to the console in the format
 * `{first,second}`. It iterates over each pair in the vector and prints the `first` and `second` elements
 * enclosed in curly braces and separated by a comma.
 *
 * @param pairs The array of pairs to be printed.
 */
void printPairArray(vector<pair<int, int>> pairs)
{
  for (auto i : pairs)
  {
    cout << "{" << i.first << "," << i.second << "} ";
  }
}

/**
 * @brief Generates a one-dimensional visualization of the secondary structure.
 *
 * This function takes an RNA sequence and an array of bonded pairs of indices, and generates a one-dimensional
 * visualization of the secondary structure of the RNA sequence. The visualization is a string where each character
 * represents a nucleotide in the sequence. Unpaired nucleotides are represented by a dot (`.`), while paired
 * nucleotides are represented by an opening parenthesis `(` for the first nucleotide in the pair, and a closing
 * parenthesis `)` for the second nucleotide in the pair.
 *
 * The function first initializes a string `oned_viz` with dots for each nucleotide in the sequence. Then, it iterates
 * over each bonded pair of indices and replaces the corresponding characters in `oned_viz` with opening and closing
 * parentheses, respectively.
 *
 * @param seq The RNA sequence as a string.
 * @param pairs The array of bonded pairs of indices.
 * @return The one-dimensional visualization string.
 */
string oneDViz(string seq, vector<pair<int, int>> pairs)
{
  string oned_viz = "";
  for (int i = 0; i < seq.size(); i++)
  {
    oned_viz.push_back('.');
  }
  for (auto i : pairs)
  {
    oned_viz[i.first] = '(';
    oned_viz[i.second] = ')';
  }
  return oned_viz;
}

/**
 * @brief Generates a comma-separated string of bonded pair indices for use with a Python script.
 *
 * This function takes an array of bonded pairs of indices and generates a comma-separated string representation
 * of these pairs. The string is formatted in a way that can be used as an argument to a Python script for further
 * processing or visualization.
 *
 * The function iterates over each pair in the input array and appends the `first` and `second` elements of the pair
 * to the output string, separated by a comma. A comma is also appended after each pair, except for the last one.
 *
 * @param pairs The array of bonded pairs of indices.
 * @return The comma-separated string of bonded pair indices.
 */
string pythonCSIndex(vector<pair<int, int>> pairs)
{
  string ans = "";
  for (auto i : pairs)
  {
    ans += to_string(i.first);
    ans.push_back(',');
    ans += to_string(i.second);
    ans.push_back(',');
  }
  ans.pop_back();
  return ans;
}

/**
 * @brief The main function of the program.
 *
 * The main function performs the following steps:
 * 1. Prompts the user to enter an RNA sequence and reads it from the input.
 * 2. Converts the input sequence to uppercase if it was entered in lowercase.
 * 3. Calls the optimalMatching function to compute the optimal secondary structure of the RNA sequence.
 * 4. Initializes a vector to store the bonded pairs of indices.
 * 5. Calls the traceback function to reconstruct the optimal secondary structure and store the bonded pairs in the vector.
 * 6. Prints the indices of the bonded pairs.
 * 7. Prints the actual bonded pairs of nucleotides.
 * 8. Calls the oneDViz function to generate a one-dimensional visualization of the secondary structure and prints it.
 * 9. Constructs a command to run a Python script for two-dimensional visualization of the secondary structure.
 * 10. Executes the Python script with the RNA sequence and the bonded pair indices as arguments.
 *
 * @return The exit status of the program.
 */
int main()
{
  string rna_seq;
  cout << "Enter RNA sequence: ";
  cin >> rna_seq;
  if (rna_seq[0] >= 'a')
  {
    for (int i = 0; i < rna_seq.size(); i++)
    {
      rna_seq[i] = (char)(rna_seq[i] - 32);
    }
  }
  vector<vector<int>> opt = optimalMatching(rna_seq);
  vector<pair<int, int>> bonded_pair_index;
  traceback(0, rna_seq.size() - 1, rna_seq, opt, bonded_pair_index);
  cout << "\nIndex of pairs: ";
  printPairArray(bonded_pair_index);
  cout << "\nThe Pairs: ";
  for (auto i : bonded_pair_index)
  {
    cout << "{" << rna_seq[i.first] << "," << rna_seq[i.second] << "} ";
  }
  string oned_viz = oneDViz(rna_seq, bonded_pair_index);
  cout << "\n\n1 Dimensional Visualization:" << endl;
  cout << oned_viz << endl;
  string pythonscript = "./twod_viz.py";
  string cmd = "python " + pythonscript;
  cmd.push_back(' ');
  cmd += rna_seq;
  cmd.push_back(' ');
  cmd += pythonCSIndex(bonded_pair_index);
  int res = system(cmd.c_str());
}