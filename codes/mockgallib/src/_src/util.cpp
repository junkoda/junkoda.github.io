
void assert_double(const char file[], const unsigned int line, const double x, const double x_expected, const double eps)
{
  double rel_error= fabs(x-x_expected)/x_expected;
  if(rel_error >= eps) {
    msg_abort("Assesion error at %s %u: %e != %e, %e error > %e required\n",
	  file, line, x, x_expected, rel_error, eps);
  }
}
