#include <string>
#include <memory>
#include <iostream>
#include <fstream>
#include <sstream>
#include <map>
#include <unordered_map>
#include <vector>
#include <stack>
#include <boost/optional.hpp>
#include <boost/program_options.hpp>
#include "four_fours.h"

using namespace std;
namespace bpo = boost::program_options;

/// 等式のキャッシュ
map<string, equations> cache;
vector<u_op> u_ops;
vector<b_op> b_ops;
map<unsigned char, op_c> op_cs;

static const double CUTOFF = 1e7;
static const double EPS = 1e-9;

/// 単項演算子の連続適用最大数
static int unary_limit = 2;
/// 平方根演算に制限を設けるか
static bool no_cutoff_sqrt = false;
/// 静かに
static bool quiet = false;

/// 演算子の準備
void setup_ops();

int main(int argc, char* argv[]) {

	string seed = "4444";
	bpo::options_description options("Options"), hidden("Hidden"), all_options("All Options");
	options.add_options()
	("help,h", "Shows this help")
	("unary-limit,u", bpo::value<int>()->default_value(2), "Maximum number of consecutive application of unary operators.")
	("quiet,q", bpo::value<bool>(), "Suppress progress reporting.")
	("no-cutoff-sqrt", "Do not restrict sqrt operand.")
	("output-file,o", bpo::value<string>(), "Output file name.");
	hidden.add_options()("numbers", bpo::value<string>());
	bpo::positional_options_description positional;
	positional.add("numbers", -1);
	all_options.add(options).add(hidden);

	bpo::variables_map vm;
	try {
		// 解析
		bpo::store(bpo::command_line_parser(argc, argv).options(all_options).positional(positional).run(), vm);
	} catch (exception &e) {
		cout << e.what() << endl;
		return 1;
	}

	// show help
	if (vm.count("help")) {
		cout << options << endl;
		return 0;
	}

	unary_limit = vm.count("unary-limit") ? vm["unary-limit"].as<int>() : 2;
	seed = vm.count("numbers") ? vm["numbers"].as<string>() : "4444";
	no_cutoff_sqrt = vm.count("no-cutoff-sqrt");

	if (seed.find_first_not_of("0123456789") != string::npos) {
		cout << "Invalid input." << endl;
		return 1;
	}

	setup_ops();
	equations&& e = think(seed);

	string output_file = vm.count("output-file") ? vm["output-file"].as<string>() : (seed + ".html");
	ofstream fout(output_file);


	fout << "<!DOCTYPE html><html><head><script type=\"text/javascript\"src=\"https://cdn.mathjax.org/mathjax/latest/MathJax.js?config=TeX-AMS-MML_HTMLorMML\"></script></head><body>" << endl;

	int c = 0;
	for (int i = 0; i <= 1000; i++) {
		auto found = e.find((val_t) i);
		if (found != e.end()) {
			c++;
			fout << "\\[" << i << "=" << exptex(found->second.rep.get(), seed) << "\\]" << endl;
		}
	}

	fout << "</body></html>" << endl;

	cout << "found: " << c << endl;
}


/// 与えられた文字列の範囲で計算します。
equations&& think(string numbers) {

	auto it = cache.find(numbers);
	if (it != cache.end()) { return move(it->second); } // もう計算しなくてよい

	if (!quiet) { cout << "-- Begin " << numbers << " --" << endl; }

	// 零項演算
	equations eqns;
	exp_t nullary;
	nullary.cost = 0;
	nullary.rep = move(n_combine(numbers.length()));
	eqns[atoi(numbers.c_str())] = move(nullary);

	for (int i = 1; i < numbers.length(); i++) {
		// 二項演算で結合（ここで並列化）
		vector<equations> b_exts;

		equations&& lhs = think(numbers.substr(0, i));
		equations&& rhs = think(numbers.substr(i));

		#pragma omp parallel for
		for (int j = 0; j < b_ops.size(); j++) {
			equations b_ext = move(b_extend(move(lhs), move(rhs), b_ops[j]));
			#pragma omp critical
			b_exts.push_back(move(b_ext));
		}
		// マージ
		for (auto it = b_exts.begin(); it != b_exts.end(); ++it) {
			for (auto jt = it->begin(); jt != it->end(); ++jt) {
				auto found = eqns.find(jt->first);
				// 見つからないか、コストが低ければ採用
				if (found == eqns.end() || found->second.cost > jt->second.cost) {
					eqns[jt->first] = move(jt->second);
				}
			}
		}
	}

	// 単項演算子で拡大
	// 既知でないものが現れたときだけ追加するようにすれば並列化できる
	vector<equations> u_steps; // 各段階の拡大
	for (int i = 0; i < unary_limit; i++) {
		vector<equations> u_op_exts;
		#pragma omp parallel for
		for (int j = 0; j < u_ops.size(); j++) {
			equations u_ext = move(u_extend(move(i == 0 ? eqns : u_steps[i - 1]), move(eqns), u_ops[j]));
			#pragma omp critical
				u_op_exts.push_back(move(u_ext));
		}
		// マージ
		equations u_eqns;
		for (auto it = u_op_exts.begin(); it != u_op_exts.end(); ++it) {
			for (auto jt = it->begin(); jt != it->end(); ++jt) {
				auto found = u_eqns.find(jt->first);
				// 見つからないか、コストが低ければ採用
				if (found == u_eqns.end() || found->second.cost > jt->second.cost) {
					u_eqns[jt->first] = move(jt->second);
				}
			}
		}

		u_steps.push_back(move(u_eqns));
	}

	// 各段階の拡大を統合
	for (auto it = u_steps.begin(); it != u_steps.end(); ++it) {
		for (auto jt = it->begin(); jt != it->end(); ++jt) {
			auto found = eqns.find(jt->first);
			// 見つからないか、コストが低ければ採用
			if (found == eqns.end() || found->second.cost > jt->second.cost) {
				eqns[jt->first] = move(jt->second);
			}
		}
	}

	if (!quiet) { cout << "-- Finish " << numbers << " --" << endl; }

	cache[numbers] = move(eqns);
	return move(cache[numbers]);
}


/// 二項演算による知識の拡大
equations b_extend(equations&& lhs, equations&& rhs, b_op op) {
	equations exts;
	for (auto itl = lhs.begin(); itl != lhs.end(); ++itl) {
		for (auto itr = rhs.begin(); itr != rhs.end(); ++itr) {
			if (auto const number = round_if_needed(op.f(itl->first, itr->first))) {
				if (!cutoff(*number)) {
					unsigned char cost = itl->second.cost + itr->second.cost + op.cost;
					auto found = exts.find(*number);
					// 見つからないか、コストが低ければ採用
					if (found == exts.end() || found->second.cost > cost) {
						exp_t e;
						e.rep = move(b_combine(itl->second.rep.get(), itr->second.rep.get(), op.token));
						e.cost = cost;
						exts[*number] = move(e);
					}
				}
			}
		}
	}
	return exts;
}

/// 単項演算による知識の拡大
/// seedから拡大し、seedにもpoolにも入っていないものを作る
equations u_extend(equations&& seed, equations&& pool, u_op op) {
	equations exts;
	for (auto it = seed.begin(); it != seed.end(); ++it) {
		if (auto const number = round_if_needed(op.f(it->first))) {
			if (cutoff(*number)) { continue; }
			unsigned char cost = it->second.cost + op.cost;
			auto seed_found = seed.find(*number);
			auto pool_found = pool.find(*number);
			if (seed_found != seed.end() && seed_found->second.cost <= cost) { continue; }
			if (pool_found != pool.end() && pool_found->second.cost <= cost) { continue; }
			exp_t e;
			e.rep = move(u_combine(it->second.rep.get(), op.token));
			e.cost = cost;
			exts[*number] = move(e);
		}
	}
	return exts;
}

inline unsigned char token_at(exp_rep_t e, int i) {
	return (i % 2 == 0) ? (e[i / 2] >> 4) : (e[i / 2] & 0x0f);
}

/// 式のトークン数を返します。
size_t explen(exp_rep_t e) {
	for (int i = 0;;i++) {
		if ((e[i] & 0x0f) == 0) {
			return (e[i] & 0xf0) == 0 ? i * 2 : i * 2 + 1;
		}
	}
}

string expstr(exp_rep_t e) {
	stringstream ss;
	for (int i = 0; i < explen(e); i++) {
		ss << (int) token_at(e, i) << ",";
	}
	return ss.str();
}

inline opt_val_t maybe_int(val_t v) {
		val_t rounded = round(v);
		return abs(rounded - v) < EPS ? opt_val_t(rounded) : opt_val_t();
}

inline opt_val_t round_if_needed(opt_val_t v) {
	if (v) {
		if (auto const i = maybe_int(*v)) {
			return i;
		} else {
			return v;
		}
	} else {
		return v;
	}
}

inline bool cutoff(val_t v) {
	return abs(v) > CUTOFF || (v != 0 && abs(v) < EPS);
}

/// 式を二項演算子で結合します。
exp_rep_p b_combine(exp_rep_t e1, exp_rep_t e2, const unsigned char t) {
	size_t l1 = explen(e1);
	size_t l2 = explen(e2);
	size_t l = l1 + l2 + 1;
	exp_rep_p e(new unsigned char[l / 2 + 1]); // lが奇数なら半分＋端数１、偶数なら半分＋ゼロ終端
	// copy of l1
	memcpy(e.get(), e1, l1 / 2 + 1);

	// copy of l2
	if (l1 % 2 == 0) {
		memcpy(e.get() + l1 / 2, e2, l2 / 2 + 1);
	} else {
		if (l2 % 2 == 0) {
			e[l1 / 2] |= (e2[0] >> 4);
			for (int i = 1; i < l2 / 2; i++) {
				e[l1 / 2 + i] = ((e2[i - 1] & 0x0f) << 4) + ((e2[i] & 0xf0) >> 4);
			}
			e[l / 2 - 1] = (e2[l2 / 2 - 1] << 4);
		} else {
			e[l1 / 2] |= (e2[0] >> 4);
			for (int i = 1; i <= l2 / 2; i++) {
				e[l1 / 2 + i] = ((e2[i - 1] & 0x0f) << 4) + ((e2[i] & 0xf0) >> 4);
			}
		}
	}

	// copy t
	if (l % 2 == 0) {
		e[l / 2 - 1] |= (t & 0x0f);
		e[l / 2] = '\0';
	} else {
		e[l / 2] = (t << 4);
	}
	return e;
}

exp_rep_p u_combine(exp_rep_t e, const unsigned char t) {
	size_t l = explen(e) + 1;
	exp_rep_p r(new unsigned char[l / 2 + 1]);
	memcpy(r.get(), e, l / 2 + 1);

	if (l % 2 == 0) {
		r[l / 2 - 1] |= (t & 0x0f);
		r[l / 2] = '\0';
	} else {
		r[l / 2] = (t << 4);
	}
	return r;
}

exp_rep_p n_combine(size_t length) {

	exp_rep_p e(new unsigned char[length / 2 + 1]);

	if (length == 0) {
		e[0] = (unsigned char) 0x00; // 終端文字のみ
	} else if (length == 1) {
		e[0] = (unsigned char) 0x10; // 数字先頭と終端文字
	} else {
		e[0] = (unsigned char) 0x12; // 数字先頭と数字継続
		for (int i = 2; i <= length - 2; i += 2) {
			e[i / 2] = (unsigned char) 0x22; // 数字継続
		}
		if (length % 2 == 0) {
			e[length / 2] = (unsigned char) 0x00; // 終端文字
		} else {
			e[length / 2] = (unsigned char) 0x20; // 数字継続と終端文字
		}
	}
	return e;
}

void setup_ops() {

	// 0x00, 0x01, 0x02は予約済み（終端文字・数字開始・数字継続）
	unsigned char token = 0x02;


	u_ops.push_back(u_op(++token, 2, [](val_t v) { return -v; }));
	op_cs[token] = op_c("-", "", 3, true);

	u_ops.push_back(u_op(++token, 4, [](val_t v) { return v >= 0 && (no_cutoff_sqrt || maybe_int(v * v * v * v)) ? opt_val_t(sqrt(v)) : opt_val_t(); }));
	op_cs[token] = op_c("\\sqrt{", "}", 1, false);

	u_ops.push_back(u_op(++token, 8, [](val_t v) { return v >= 0 && (no_cutoff_sqrt || maybe_int(v * v)) ? opt_val_t(sqrt(sqrt(v))) : opt_val_t(); }));
	op_cs[token] = op_c("\\sqrt{\\sqrt{", "}}", 1, false);

	u_ops.push_back(u_op(++token, 6, [](val_t v) {
		static int table[11] = {1, 1, 2, 6, 24, 120, 720, 5040, 40320, 362880, 3628800};
		if (auto const integer = maybe_int(v)) {
			if (*integer >= 0 && *integer <= 10) {
				return opt_val_t((val_t) table[(int) *integer]);
			}
		}
		return opt_val_t();
	}));
	op_cs[token] = op_c("", "!", 2, true);


	b_ops.push_back(b_op(++token, 1, [](val_t lhs, val_t rhs) { return lhs + rhs; }));
	op_cs[token] = op_c("", "+", "", 6, true, true, true, true);

	b_ops.push_back(b_op(++token, 2, [](val_t lhs, val_t rhs) { return lhs - rhs; }));
	op_cs[token] = op_c("", "-", "", 6, true, false, true, true);

	b_ops.push_back(b_op(++token, 3, [](val_t lhs, val_t rhs) { return lhs * rhs; }));
	op_cs[token] = op_c("", "\\times", "", 5, true, true, true, true);

	b_ops.push_back(b_op(++token, 4, [](val_t lhs, val_t rhs) { return rhs != 0 ? opt_val_t(lhs / rhs) : opt_val_t(); }));
	op_cs[token] = op_c("\\frac{", "}{", "}", 5, true, false, false, false);

	b_ops.push_back(b_op(++token, 6, [](val_t lhs, val_t rhs) { return maybe_int(rhs) && lhs >= 0 ? opt_val_t(pow(lhs, rhs)) : opt_val_t(); }));
	op_cs[token] = op_c("{", "}^{", "}", 3, false, true, true, false);
}

string exptex(exp_rep_t e, string seed) {
	stack<pair<string, int> > st; // （部分式文字列、最後に追加された演算子の優先度）のペアのスタック
	string paren_left = "\\left("; string paren_right = "\\right)";
	auto sit = seed.begin();
	for (int i = 0; i < explen(e); i++) {
		unsigned char token = token_at(e, i);
		if (token == 0x01) {
			st.push(make_pair(string(1, *sit++), 0));
		} else if (token == 0x02) {
			string subexp = st.top().first + string(1, *sit++);
			st.pop();
			st.push(make_pair(subexp, 0));
		} else {
			auto found = op_cs.find(token);
			if (found == op_cs.end()) { continue; }
			auto op_c = found->second;

			if (op_c.unary) {
				auto opr = st.top();
				string subexp = opr.first;
				st.pop();
				if (op_c.needs_paren && opr.second >= op_c.precedence) { subexp = paren_left + subexp + paren_right; }
				subexp = op_c.prefix + subexp + op_c.suffix;
				st.push(make_pair(subexp, op_c.precedence));
			} else {
				auto op2 = st.top();
				st.pop();
				auto op1 = st.top();
				st.pop();
				string subexp1 = op1.first;
				string subexp2 = op2.first;

				if (op_c.needs_paren_op1 && (op1.second > op_c.precedence || (op1.second == op_c.precedence && !op_c.left_associative))) {
					subexp1 = paren_left + subexp1 + paren_right;
				}
				if (op_c.needs_paren_op2 && (op2.second > op_c.precedence || (op2.second == op_c.precedence && !op_c.right_associative))) {
					subexp2 = paren_left + subexp2 + paren_right;
				}

				string subexp = op_c.prefix + subexp1 + op_c.infix + subexp2 + op_c.suffix;
				st.push(make_pair(subexp, op_c.precedence));
			}
		}
	}

	return st.top().first;
}
