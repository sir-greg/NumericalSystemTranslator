#include "mainwindow.h"
#include "./ui_mainwindow.h"
#include <iostream>
#include <string>
#include <vector>
#include <algorithm>
#include <set>

using std::vector;
using std::to_string;
using std::string;
using std::cout;
using std::cin;
using std::set;

namespace LongArithmetic {

    vector<char> digs {
        '0', '1', '2', '3', '4', '5', '6', '7', '8', '9',
        'A', 'B', 'C', 'D', 'E', 'F', 'G', 'H', 'I', 'J',
        'K', 'L', 'M', 'N', 'O', 'P', 'Q', 'R', 'S', 'T',
        'U', 'V', 'W', 'X', 'Y', 'Z'
    };

    class bigNumInt {

    private:
        vector<int64_t> digits = {};
        static const int64_t dig_per_ind = 9;
        int64_t real_base = 10;
        static const int64_t base = 1000000000LL;
        int sign = 1;

        int64_t div(int64_t a, int64_t b) {
            if (a > 0) {
                return a / b;
            } else {
                int sign = 0;
                if (b > 0) {
                    sign = -1;
                } else {
                    sign = 1;
                }
                a = abs(a);
                b = abs(b);
                return (a / b + (a % b != 0)) * sign;
            }
        }

        int64_t mod(int64_t a, int64_t b) {
            return (a % b + b) % b;
        }

    public:
        bigNumInt() {}

        explicit bigNumInt(const string& str, int64_t str_base = 10) {
            int64_t part;
            real_base = str_base;
            if (str[0] == '-') {
                sign = -1;
                part = 1;
            } else {
                sign = 1;
                part = 0;
            }
            if (str_base == 10) {
                digits.resize((str.size() + dig_per_ind - 1) / dig_per_ind);
                for (int64_t ind = str.size(); ind > part; ind -= dig_per_ind) {
                    if (ind < dig_per_ind) {
                        digits[ind / dig_per_ind - (ind % dig_per_ind == 0)] =
                            stoll(str.substr(part, ind - part));
                    } else {
                        digits[ind / dig_per_ind - (ind % dig_per_ind == 0)] =
                            stoll(str.substr(ind - dig_per_ind, dig_per_ind));
                    }
                }
            } else {
                digits = { 0 };
                for (char c : str) {
                    if (c == '-') continue;
                    *this *= str_base;
                    if (c <= '9') {
                        *this += bigNumInt(to_string(c - '0'));
                    } else {
                        *this += bigNumInt(to_string(c - 'A' + 10));
                    }
                }
            }
        }

        explicit bigNumInt(const int64_t number) {
            *this = bigNumInt(to_string(number));
        }

        bigNumInt(const bigNumInt& other) {
            this->digits = other.digits;
            this->sign = other.sign;
        }

        int64_t get_digit(const int index) const {
            return digits[index];
        }

        int64_t& operator[](const int index) {
            return digits[index];
        }

        int64_t get_base() const {
            return real_base;
        }

        int64_t make_int64_t() {
            int64_t res = this->get_digit(0);
            if (this->get_size() > 1) {
                res *= base;
                res += this->get_digit(1);
            }
            return res;
        }

        int get_sign() const {
            return sign;
        }

        void set_sign(int new_sign) {
            sign = new_sign;
        }

        void set_size(int64_t size) {
            int64_t times = size - get_size();
            while (times--) {
                digits.insert(digits.begin(), 0);
            }
        }

        int64_t get_size() const {
            return digits.size();
        }

        void delete_leading_zeros() {
            while (get_size() > 1 && *digits.begin() == 0) {
                digits.erase(digits.begin());
            }
        }

        void add_zeros(int64_t zeros) {
            while (zeros--) {
                digits.push_back(0);
            }
        }

        vector<int64_t>::iterator end() {
            return digits.end();
        }

        vector<int64_t>::iterator begin() {
            return digits.begin();
        }

        bigNumInt operator-() const {
            bigNumInt res(*this);
            res.switch_sign();
            return res;
        }

        void switch_sign() {
            this->sign *= -1;
        }

        bool operator>=(bigNumInt other) const {
            if (get_sign() > other.get_sign()) {
                return true;
            } else if (get_sign() < other.get_sign()) {
                return false;
            } else if (get_sign() == -1) {
                return !((-*this) >= (-other));
            } else if (get_size() == other.get_size()) {
                for (int64_t ind = 0; ind < get_size(); ++ind) {
                    if (get_digit(ind) != other.get_digit(ind)) {
                        return get_digit(ind) > other.get_digit(ind);
                    }
                }
                return true;
            }
            return get_size() > other.get_size();
        }

        bool operator>(bigNumInt other) const {
            if (get_sign() > other.get_sign()) {
                return true;
            } else if (get_sign() < other.get_sign()) {
                return false;
            } else if (get_sign() == -1) {
                return !((-*this) > (-other));
            } else if (get_size() == other.get_size()) {
                for (int64_t ind = 0; ind < get_size(); ++ind) {
                    if (get_digit(ind) != other.get_digit(ind)) {
                        return get_digit(ind) > other.get_digit(ind);
                    }
                }
                return false;
            }
            return get_size() > other.get_size();
        }

        bool operator<(bigNumInt other) const {
            return !(*this >= other);
        }

        bool operator<=(bigNumInt other) const {
            return !(*this > other);
        }

        bool operator==(bigNumInt other) const {
            return !(*this > other) && !(*this < other);
        }

        bool operator!=(bigNumInt other) const {
            return !(*this == other);
        }

        bigNumInt operator=(const bigNumInt& other) {
            this->digits = other.digits;
            this->sign = other.sign;
            return *this;
        }

        bigNumInt operator+(bigNumInt other) {
            bool switch_sign = 0;
            if (this->get_sign() == -1) {
                if (other.get_sign() == -1) {
                    switch_sign = 1;
                    this->switch_sign();
                    other.switch_sign();
                } else {
                    if (-*this > other) {
                        switch_sign = 1;
                        this->switch_sign();
                        other.switch_sign();
                    }
                }
            } else {
                if (other.get_sign() == -1 && -other > *this) {
                    switch_sign = 1;
                    this->switch_sign();
                    other.switch_sign();
                }
            }
            int64_t common_size = std::max(this->get_size(), other.get_size()) + 1;
            this->set_size(common_size);
            other.set_size(common_size);
            bigNumInt res;
            res.set_size(common_size);
            int64_t carry = 0;
            int64_t dig = 0;
            for (int64_t ind = common_size - 1; ind > -1; --ind) {
                dig = this->get_digit(ind) * this->get_sign() +
                    other.get_digit(ind) * other.get_sign() + carry;
                res[ind] = mod(dig, base);
                carry = div(dig, base);
            }
            res.delete_leading_zeros();
            this->delete_leading_zeros();
            if (switch_sign) {
                res.switch_sign();
                this->switch_sign();
            }
            return res;
        }

        bigNumInt operator+=(bigNumInt other) {
            *this = *this + other;
            return *this;
        }

        bigNumInt operator-(bigNumInt other) {
            return *this + (-other);
        }

        bigNumInt operator-=(bigNumInt other) {
            *this = *this - other;
            return *this;
        }

        void operator++() {
            *this += bigNumInt("1");
        }

        void operator--() {
            *this -= bigNumInt("1");
        }

        bigNumInt operator++(int) {
            *this += bigNumInt("1");
            return *this - bigNumInt("1");
        }

        bigNumInt operator--(int) {
            *this -= bigNumInt("1");
            return *this + bigNumInt("1");
        }

        bigNumInt operator*(int64_t number) {
            if (!this->count_digits()) return bigNumInt(0);
            bigNumInt res;
            res.set_size(this->get_size() + 1);
            bool switch_sign = 0;
            if (number < 0) {
                number = -number;
                switch_sign = 1;
            }
            bool this_switched = 0;
            if (this->get_sign() == -1) {
                this->switch_sign();
                this_switched = 1;
                switch_sign = !switch_sign;
            }
            int64_t dig;
            int64_t carry = 0;
            for (int64_t ind = res.get_size() - 1; ind > 0; --ind) {
                dig = this->get_digit(ind - 1) * number + carry;
                res[ind] = mod(dig, base);
                carry = div(dig, base);
            }
            if (carry) {
                res[0] = carry;
            }
            res.delete_leading_zeros();
            if (switch_sign) {
                res.switch_sign();
            }
            if (this_switched) {
                this->switch_sign();
            }
            return res;
        }

        bigNumInt operator*=(int64_t number) {
            *this = *this * number;
            return *this;
        }

        bigNumInt operator*(bigNumInt other) {
            if (this->count_digits() * other.count_digits() == 0) return bigNumInt(0);
            bigNumInt res, extraNum;
            bool switch_sign = 0;
            bool this_switched = 0;
            if (other.get_sign() == -1) {
                other.switch_sign();
                switch_sign = !switch_sign;
            }
            if (this->get_sign() == -1) {
                this->switch_sign();
                this_switched = 1;
                switch_sign = !switch_sign;
            }
            for (int64_t ind = this->get_size() - 1; ind > -1; --ind) {
                extraNum = other * this->get_digit(ind);
                extraNum.add_zeros(this->get_size() - 1 - ind);
                res += extraNum;
            }
            if (switch_sign) {
                res.switch_sign();
            }
            if (this_switched) {
                this->switch_sign();
            }
            return res;
        }

        bigNumInt operator*=(bigNumInt other) {
            *this = *this * other;
            return *this;
        }

        bigNumInt operator/(int64_t number) {
            bigNumInt res;
            int64_t carry = 0;
            int64_t dig = 0;
            bool switch_sign = 0;
            bool this_switched = 0;
            if (number < 0) {
                number = -number;
                switch_sign = 1;
            }
            if (this->get_sign() == -1) {
                this->switch_sign();
                this_switched = 1;
                switch_sign = !switch_sign;
            }
            for (int64_t ind = 0; ind < this->get_size(); ++ind) {
                res.add_zeros(1);
                res[ind] = (carry + this->get_digit(ind)) / number;
                carry = (carry + this->get_digit(ind)) % number * base;
            }
            res.delete_leading_zeros();
            if (switch_sign) {
                res.switch_sign();
                number = -number;
            }
            if (this_switched) {
                this->switch_sign();
            }
            if (res.get_sign() == -1 && *this % number != 0) {
                --res;
            }
            return res;
        }

        bigNumInt operator/=(int64_t number) {
            *this = *this / number;
            return *this;
        }

        int64_t operator%(int64_t number) {
            bigNumInt res;
            int64_t carry = 0;
            int64_t dig = 0;
            bool switch_sign = 0;
            bool this_switched = 0;
            if (number < 0) {
                number = -number;
                switch_sign = !switch_sign;
            }
            if (this->get_sign() == -1) {
                this->switch_sign();
                switch_sign = !switch_sign;
                this_switched = 1;
            }
            for (int64_t ind = 0; ind < this->get_size(); ++ind) {
                res.add_zeros(1);
                res[ind] = (carry + this->get_digit(ind)) / number;
                carry = (carry + this->get_digit(ind)) % number * base;
            }
            carry /= base;
            if (switch_sign) {
                carry = (number - carry) % number;
            }
            if (this_switched) {
                this->switch_sign();
            }
            return carry;
        }

        bigNumInt operator/(bigNumInt other) {
            bigNumInt res, piece;
            bool switch_sign = 0;
            bool this_switch = 0;
            if (other.get_sign() == -1) {
                other.switch_sign();
                switch_sign = 1;
            }
            if (this->get_sign() == -1) {
                this->switch_sign();
                this_switch = 1;
                switch_sign = !switch_sign;
            }
            for (int64_t ind = 0; ind < this->get_size(); ++ind) {
                piece = piece * base + bigNumInt(to_string(this->get_digit(ind)));
                if (piece >= other) {
                    int64_t left = 0, right = base, mid;
                    while (right - left > 1) {
                        mid = left + (right - left) / 2;
                        if (other * mid > piece) {
                            right = mid;
                        } else {
                            left = mid;
                        }
                    }
                    res.add_zeros(1);
                    res[ind] = left;
                    piece = piece - other * left;
                } else {
                    res.add_zeros(1);
                    res[ind] = 0;
                }
            }
            res.delete_leading_zeros();
            if (switch_sign) {
                res.switch_sign();
            }
            if (this_switch) {
                this->switch_sign();
            }
            if (res.get_sign() == -1 && piece != bigNumInt("0")) {
                --res;
            }
            return res;
        }

        bigNumInt mod(bigNumInt other) {
            bigNumInt remainder = *this - (*this / other) * other;
            return remainder;
        }

        bigNumInt operator%(bigNumInt other) {
            bigNumInt remainder = this->mod(other);
            if (remainder.get_sign() == -1) {
                remainder += other;
            }
            return remainder;
        }

        bigNumInt operator/=(bigNumInt other) {
            *this = *this / other;
            return *this;
        }

        bigNumInt operator%=(bigNumInt other) {
            *this = *this % other;
            return *this;
        }

        bigNumInt sqr() {
            return *this * *this;
        }

        void change_base(const int64_t newBase) {
            this->real_base = newBase;
        }

        int64_t count_digits() {
            return (this->get_size()) ? to_string(this->get_digit(0)).size()
                + (this->get_size() - 1) * 9 : 0;
        }

        bool isEmpty() const {
            return this->get_size() == 0;
        }

        string get_str() const {
            string out;
            if (this->get_size() == 0) return out;
            if (this->get_base() == 10) {
                if (this->get_sign() == -1) {
                    out += "-";
                }
                out += to_string(this->get_digit(0));
                for (int64_t ind = 1; ind < this->digits.size(); ++ind) {
                    string output = to_string(this->digits[ind]);
                    for (int zeros = 0; zeros < dig_per_ind - output.size();
                        ++zeros) {
                        out += "0";
                    }
                    out += output;
                }
            }
            else {
                if (this->get_sign() == -1) {
                    out += "-";
                }
                string output;
                bigNumInt reserve = *this;
                while (reserve > bigNumInt(0)) {
                    output += digs[reserve % this->get_base()];
                    reserve /= this->get_base();
                }
                reverse(output.begin(), output.end());
                if (output.size() == 0) output = "0";
                out += output;
            }
            return out;
        }

        friend std::ostream& operator<<(std::ostream& os, const bigNumInt& num) {
            if (num.get_size() == 0) return os;
            if (num.get_base() == 10) {
                if (num.get_sign() == -1) {
                    os << "-";
                }
                os << num.get_digit(0);
                for (int64_t ind = 1; ind < num.digits.size(); ++ind) {
                    string output = to_string(num.digits[ind]);
                    for (int zeros = 0; zeros < dig_per_ind - output.size();
                        ++zeros) {
                        os << 0;
                    }
                    os << output;
                }
            } else {
                if (num.get_sign() == -1) {
                    os << "-";
                }
                string output;
                bigNumInt reserve = num;
                while (reserve > bigNumInt(0)) {
                    output += digs[reserve % num.get_base()];
                    reserve /= num.get_base();
                }
                reverse(output.begin(), output.end());
                if (output.size() == 0) output = "0";
                os << output;
            }
            return os;
        }

        friend std::istream& operator>>(std::istream& is, bigNumInt& num) {
            string number;
            is >> number;
            num = bigNumInt(number);
            return is;
        }
    };
    void Swap(bigNumInt& a, bigNumInt& b);
    void Swap(bigNumInt& a, bigNumInt& b) {
        bigNumInt t;
        t = a;
        a = b;
        b = t;
    }

    void Abs(bigNumInt& a);
    void Abs(bigNumInt& a) {
        if (a.get_sign() == -1) {
            a.switch_sign();
        }
    }

    bigNumInt Gcd(bigNumInt a, bigNumInt b);
    bigNumInt Gcd(bigNumInt a, bigNumInt b) {
        if (a.get_sign() == -1) {
            a.switch_sign();
        }
        if (b.get_sign() == -1) {
            b.switch_sign();
        }
        while (b != bigNumInt(0)) {
            a = a % b;
            Swap(a, b);
        }
        return a;
    }

    bigNumInt Pow(bigNumInt number, int64_t pw);
    bigNumInt Pow(bigNumInt number, int64_t pw) {
        if (pw == 0) {
            return bigNumInt(1);
        } else if (pw % 2 == 0) {
            bigNumInt res = Pow(number, pw / 2);
            return res * res;
        } else {
            bigNumInt res = Pow(number, pw / 2);
            return res * res * number;
        }
    }

    class fraction {
    private:
        bigNumInt numerator;
        bigNumInt denominator;
        int64_t base = 10;
        int sign = 1;

    public:
        fraction() {}

        fraction(bigNumInt num, bigNumInt denom, int sign_) : numerator(num), denominator(denom), sign(sign_) {}

        fraction(const string& str, const int64_t str_base = 10) {
            if (str.empty()) {
                return;
            }
            bigNumInt intPart, preperiod, period;
            int64_t intPartDigits = 0, preperiodDigits = 0, periodDigits = 0;
            int64_t point;
            int64_t start;
            int64_t zeros;
            if (str[0] == '-') {
                start = 1;
                point = 1;
                sign = -1;
            } else {
                start = 0;
                point = 0;
                sign = 1;
            }
            base = str_base;

            // reading int part

            while (point < str.size() && str[point] != '.') {
                ++point;
                ++intPartDigits;
            }

            intPart = bigNumInt(str.substr(start, point - start), str_base);
            ++point;
            start = point;

            // reading preperiod (if there is one)

            zeros = 0;

            while (point < str.size() && str[point] == '0') {
                ++zeros;
                ++point;
            }

            while (point < str.size() && str[point] != '(') {
                ++point;
            }

            if (point == start) {
                preperiod = bigNumInt();
                preperiodDigits = 0;
            } else {
                preperiod = bigNumInt(str.substr(start, point - start), str_base);
                preperiodDigits = point - start;
            }
            ++point;
            start = point;

            // reading period (if there is one)

            zeros = 0;

            while (point < str.size() && str[point] == '0') {
                ++zeros;
                ++point;
            }

            while (point < str.size() && str[point] != ')') {
                ++point;
            }

            if (point == start) {
                period = bigNumInt(0);
                periodDigits = zeros;
            } else {
                period = bigNumInt(str.substr(start, point - start), str_base);
                periodDigits = point - start;
            }

            // creating fraction
            if (periodDigits > 0) {
                numerator =
                    intPart * Pow(bigNumInt(str_base), preperiodDigits) * (Pow(bigNumInt(str_base), periodDigits) - bigNumInt(1)) +
                    preperiod * (Pow(bigNumInt(str_base), (periodDigits)) - bigNumInt(1)) + period;

                denominator =
                    Pow(bigNumInt(str_base), preperiodDigits) * (Pow(bigNumInt(str_base), periodDigits) - bigNumInt(1));
            } else {
                numerator =
                    intPart * Pow(bigNumInt(str_base), preperiodDigits) + preperiod;

                denominator =
                    Pow(bigNumInt(str_base), preperiodDigits);
            }

            bigNumInt num_denom_gcd = Gcd(numerator, denominator);

            numerator /= num_denom_gcd;

            denominator /= num_denom_gcd;

            numerator.change_base(str_base);

            denominator.change_base(str_base);
        }

        int64_t get_base() const  {
            return base;
        }

        bigNumInt get_numerator() const  {
            return numerator;
        }

        bigNumInt get_denominator() const  {
            return denominator;
        }

        int get_sign() const {
            return sign;
        }

        string get_str() const {
            string out;
            if (this->get_sign() == -1) {
                out += "-";
            }
            bigNumInt intPart = this->get_numerator() / this->get_denominator();
            intPart.change_base(this->get_base());
            out += intPart.get_str();
            bigNumInt num_left = this->get_numerator() % this->get_denominator();
            if (num_left > bigNumInt(0)) {
                int64_t oldSize = -1;
                set<bigNumInt> fractions;
                vector<bigNumInt> fractions_enum;
                string fraction_out;

                fractions.insert(num_left);
                fractions_enum.push_back(num_left);

                while (oldSize != fractions.size() && num_left != bigNumInt(0)) {
                    oldSize = fractions.size();
                    num_left *= this->get_base();
                    fraction_out += digs[(num_left / this->get_denominator()).make_int64_t()];
                    num_left %= this->get_denominator();
                    fractions.insert(num_left);
                    fractions_enum.push_back(num_left);
                }
                int64_t preperiod_endpos = 0;
                while (preperiod_endpos < fractions_enum.size() &&
                    fractions_enum[preperiod_endpos] != fractions_enum.back()) {
                    ++preperiod_endpos;
                }
                if (preperiod_endpos == fractions_enum.size() - 1) {
                    if (!fraction_out.empty()) {
                        out += ".";
                        out += fraction_out;
                    }
                }
                else {
                    out += "." + fraction_out.substr(0, preperiod_endpos);
                    string period_str = fraction_out.substr(preperiod_endpos);
                    // period_str.pop_back();
                    out += "(" + period_str + ")";
                }
            }
            return out;
        }

        void change_base(const int64_t newBase) {
            base = newBase;
            numerator.change_base(base);
            denominator.change_base(base);
        }

        friend std::ostream& operator<<(std::ostream& os, const fraction& Fraction) {
            if (Fraction.get_sign() == -1) {
                os << "-";
            }
            bigNumInt intPart = Fraction.get_numerator() / Fraction.get_denominator();
            intPart.change_base(Fraction.get_base());
            bigNumInt num_left = Fraction.get_numerator() % Fraction.get_denominator();
            if (num_left > bigNumInt(0)) {
                int64_t oldSize = -1;
                set<std::pair<bigNumInt, bigNumInt>> fractions;
                vector<std::pair<bigNumInt, bigNumInt>> fractions_enum;
                string fraction_out;

                fractions.insert({ num_left, Fraction.get_denominator() });
                fractions_enum.push_back({ num_left, Fraction.get_denominator() });

                while (oldSize != fractions.size() && num_left != bigNumInt(0)) {
                    oldSize = fractions.size();
                    num_left *= Fraction.get_base();
                    fraction_out += digs[(num_left / Fraction.get_denominator()).make_int64_t()];
                    num_left %= Fraction.get_denominator();
                    fractions.insert({ num_left, Fraction.get_denominator() });
                    fractions_enum.push_back({ num_left, Fraction.get_denominator() });
                }
                int64_t preperiod_endpos = 0;
                while (preperiod_endpos < fractions_enum.size() &&
                    fractions_enum[preperiod_endpos] != fractions_enum.back()) {
                    ++preperiod_endpos;
                }
                if (preperiod_endpos == fractions_enum.size() - 1) {
                    if (!fraction_out.empty()) {
                        os << ".";
                        os << fraction_out;
                    }
                }
                else {
                    os << "." << fraction_out.substr(0, preperiod_endpos);
                    string period_str = fraction_out.substr(preperiod_endpos);
                    // period_str.pop_back();
                    os << "(" << period_str << ")";
                }
            }
            return os;
        }

        friend std::istream& operator>>(std::istream& is, fraction& Fraction) {
            string str;
            int64_t base;
            is >> str >> base;
            Fraction = fraction(str, base);
            return is;
        }
    };

}

MainWindow::MainWindow(QWidget *parent)
    : QMainWindow(parent)
    , ui(new Ui::MainWindow)
{
    ui->setupUi(this);
    QObject::connect(ui->CalculateTheResultButton, &QPushButton::clicked,
                     this, &MainWindow::CalculateTheResult);
}

MainWindow::~MainWindow()
{
    delete ui;
}

void MainWindow::CalculateTheResult()
{
    using LongArithmetic::fraction;
    bool flag;
    QString number_to_translate = ui->NumberLineEdit->text();
    qint32 from_base = ui->FromBaseLineEdit->text().toInt(&flag);
    if (!flag) {
        return;
    }
    qint32 to_base = ui->ToBaseLineEdit->text().toInt(&flag);
    if (!flag) {
        return;
    }

    fraction fraction_number(number_to_translate.toStdString(), from_base);
    fraction_number.change_base(to_base);
    ui->ResultTextEdit->setText(QString::fromStdString(fraction_number.get_str()));

}

