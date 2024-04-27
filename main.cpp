#include <bits/stdc++.h>
using namespace std;

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

vector<vector<int>> optimalMatching(string str)
{
  int n = str.size();
  // map<pair<int,int>,vector<pair<int,int>> max_cuzof;
  vector<vector<int>> opt(n, vector<int>(n, -1));
  for (int diff = 0; diff <= min(4, n - 1); diff++)
  {
    for (int i = 0; i < n - diff; i++)
    {
      int j = i + diff;
      opt[i][j] = 0;
      opt[j][i] = 0;
    }
  }
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
      // cout<<opt[i][j]<<" "<<i<<","<<j<<endl;
    }
  }
  // printDP(opt);
  return opt;
}

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

void printPairArray(vector<pair<int, int>> pairs)
{
  for (auto i : pairs)
  {
    cout << "{" << i.first << "," << i.second << "} ";
  }
}

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