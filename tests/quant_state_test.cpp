// Quantifier state across repeated and mixed initialisation. Build with -fsanitize=address so that leaks and
// double frees fail the run too (see run_quant_tests.sh).
// usage: quant_state_test <index prefix> <protein index prefix> <classification file> <nodes.dmp> <names.dmp>
#include <stdio.h>
#include <string.h>
#include <string>
#include "Quantifier.hpp"

static std::string render(Quantifier &q, int format = 0, bool emReads = false)
{
  char *buf = NULL ;
  size_t len = 0 ;
  FILE *fp = open_memstream(&buf, &len) ;
  q.Output(fp, format, emReads) ;
  fclose(fp) ;
  std::string s(buf, len) ;
  free(buf) ;
  return s ;
}

static int failures = 0 ;
static void check(const char *name, bool ok, const std::string &a = "", const std::string &b = "")
{
  printf("%-66s %s\n", name, ok ? "ok" : "FAIL") ;
  if (!ok)
  {
    printf("--- expected\n%s--- got\n%s", a.c_str(), b.c_str()) ;
    ++failures ;
  }
}

int main(int argc, char **argv)
{
  if (argc < 6)
  {
    fprintf(stderr, "usage: %s <index prefix> <protein index prefix> <classification> <nodes.dmp> <names.dmp>\n", argv[0]) ;
    return 2 ;
  }
  char *idx = argv[1], *protIdx = argv[2], *cls = argv[3], *nodes = argv[4], *names = argv[5] ;
  setvbuf(stdout, NULL, _IONBF, 0) ;   // keep the results printed before a crash

  Quantifier once ;
  once.Init(idx) ;
  once.LoadReadAssignments(cls, 0, 0, 0) ;
  once.Quantification() ;
  std::string a = render(once) ;

  Quantifier twice ;
  twice.Init(idx) ;
  twice.Init(idx) ;
  twice.LoadReadAssignments(cls, 0, 0, 0) ;
  twice.Quantification() ;
  twice.Quantification() ;
  std::string b = render(twice) ;
  check("index: Init x2 + Quantification x2 equals one of each", a == b, a, b) ;

  Quantifier dump ;
  dump.Init(nodes, names, NULL) ;
  dump.LoadReadAssignments(cls, 0, 0, 0) ;
  dump.Quantification() ;
  std::string c = render(dump) ;

  std::string d ;
  {                                   // own scope: its destructor (where a double free would show) runs here
    Quantifier mixed ;
    mixed.Init(idx) ;
    mixed.LoadReadAssignments(cls, 0, 0, 0) ;
    mixed.Init(nodes, names, NULL) ;
    mixed.LoadReadAssignments(cls, 0, 0, 0) ;
    mixed.Quantification() ;
    d = render(mixed) ;
  }
  check("index Init then taxonomy-file Init equals taxonomy-file Init alone", c == d, c, d) ;

  Quantifier stale ;
  stale.Init(nodes, names, NULL) ;
  stale.LoadReadAssignments(cls, 0, 0, 0) ;
  stale.Init(nodes, names, NULL) ;
  stale.Quantification() ;            // the assignments loaded before Init must be gone
  std::string e = render(stale) ;
  check("re-Init drops assignments loaded before it (header only)", e.find('\n') == e.size() - 1, "header line only\n", e) ;

  Quantifier none ;
  none.Quantification() ;             // must report and return, not crash
  std::string f = render(none) ;
  check("Quantification and Output before Init write nothing", f.empty(), "", f) ;

  Quantifier prot ;
  prot.Init(idx) ;
  bool ntIndex = prot.IsProtein() ;
  prot.Init(protIdx) ;
  bool aaIndex = prot.IsProtein() ;
  prot.Init(nodes, names, NULL) ;
  bool afterFiles = prot.IsProtein() ;
  prot.SetProtein(true) ;
  bool afterSet = prot.IsProtein() ;
  check("protein: read from .4.cfr, reset by taxonomy-file Init, SetProtein",
      !ntIndex && aaIndex && !afterFiles && afterSet, "0 1 0 1\n",
      std::to_string(ntIndex) + " " + std::to_string(aaIndex) + " " + std::to_string(afterFiles) + " " + std::to_string(afterSet) + "\n") ;

  return failures ? 1 : 0 ;
}
