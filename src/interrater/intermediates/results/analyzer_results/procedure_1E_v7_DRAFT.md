# 1E — Results: procedure

## 1. Purpose

Write the vaccination results of each included study into `data_extraction.md`, each with the
page it comes from. Then generate the CSV, `study_results.csv`, from `data_extraction.md`.

## 2. What this stage starts from

- `data_extraction.md` as sealed in stage 1D (`data_extraction_2026-11-09.md`).
- The study tables sealed in stage 0 (`V2_P0a_Supplement_2026-10-21.docx`), which cover 37 of
  the 113 included studies.
- The PDFs of the 113 included studies.

## 3. Who does what

- **Extraction:** a separate Claude agent for each study. Each has seen none of the earlier
  work and no other study. Claude Opus 5.5 is preferred; the model that ran is recorded.
  The agent works in two steps:
  - **Step A.** It is given this procedure and the PDF, and nothing else. It writes the
    results. They are saved as written.
  - **Step B.** It is then shown the study's existing `Outcome:` line in `data_extraction.md`.
    Where that line and the new results differ, it goes back to the paper and decides which
    is right.
- **Page check:** a script, during the extraction (section 6).
- **Consistency check:** a second such agent, after the extraction is finished, reads the
  new results and checks that they follow this procedure (section 7).
- **Comparison with the stage 0 tables:** another such agent, after the extraction is
  finished (section 7).
- **Blind verification:** the author, and a Claude Fable 5.2 agent that has seen none of
  the earlier work (section 8).
- **Decisions:** anything an agent cannot settle from the paper is held for the author.

## 4. Which results to extract

1. **Only vaccination results, of three kinds:** at least a stated number of doses (most
   often one); completion of the series; completion on time. Not extracted: predictors, odds
   ratios, costs, vaccine effectiveness, missed opportunities.
2. **The result for the whole study population**, as in the abstract or the main results
   table. If the paper gives no result for its whole population, and reports only by
   groups (age groups, periods, plans):
   - where it prints the size of every group, and each group's percentage is the plain
     share of that group's people, the result for the whole population is calculated from
     them and marked "our calculation" (section 5);
   - where a group's figure is an average or an estimate (an average of monthly figures, a
     Kaplan-Meier estimate), nothing is calculated: take the group with the largest number
     of people and say so (section 5);
   - where the sizes of the groups could not be found, take the group the paper puts
     forward as its main finding in its abstract or conclusion; if it names none, the most
     recent; and say so (section 5).
   The groups themselves are subgroups (rule 6). Where the paper gives different figures
   for the same result:
   first, take the figure the paper itself states as the result, in its abstract or text;
   if it states none, take the one with more information, for example the one whose
   numerator and denominator are clear; if they are equal, either.
3. **A separate result where the paper reports separately by:**
   - the denominator: everyone eligible, people who received a first dose, or people who
     completed the series;
   - the payer: commercial, Medicaid, Medicare (Medicare Advantage and fee-for-service
     separately where the paper separates them), other;
   - the data source, where the paper reports more than one claims database: one result
     for each, and one for all sources together where the paper states one, even without
     its numbers. Databases are never added together by us, because the same person can
     be in more than one;
   - the vaccine, where a study covers several, including a vaccine reported only for
     comparison;
   - the brand of the vaccine: one result for each brand, and one for all brands together
     where the paper gives it;
   - the schedule, where the paper judges the series against more than one (for example
     the package insert and the ACIP schedule);
   - cohorts with different observation requirements, including cohorts the paper sets up
     to test the requirement.
4. **One result, one figure:** one group of people, one outcome, one time limit, one
   percentage, with its numerator and denominator where printed. Each result becomes one
   row of the CSV. Where the paper prints only the opposite (for example the number not
   vaccinated), the result is worked out from the printed numbers and marked "our
   calculation". Where the paper prints a percentage without its numerator, and prints the
   counts that make the numerator up, the numerator is their sum, marked "our calculation",
   provided the sum gives the percentage the paper prints. If it does not, write "numerator
   not found" and hold it for the author.
5. **Several time limits.** Where the same people are read at several ages or dates, take
   the one the paper treats as its main result; if it names none, the longest. The others
   are listed, without figures, under "Other time limits reported". Where the paper
   requires longer enrolment for the later age or date, so that the number of people
   differs, each is a separate cohort and a separate result (rule 3). A different schedule
   is not another time limit; it is a separate result (rule 3).
6. **Subgroups** (age, sex, race, region, year and the like) are not extracted. Where the
   paper gives one of the three kinds of result separately for subgroups, those subgroups
   are listed, without figures, under "Subgroups explored". Characteristics the paper only
   describes, without a vaccination result for each, are not listed.
7. **Tests of the observation requirement.** Where a paper tests its observation
   requirement, each cohort is a result (rule 3). The comparison between them (a
   difference, a ratio, a sensitivity bound) is not a result. It stays in
   `Requirement effect tested`.
8. **A study with no such result** gets one line saying so and what it reports instead.
   Results are taken from the PDF only. A result the paper says is in a supplement that is
   not in the PDF is not written; it is held for the author in the words "in a supplement
   not in the PDF:", with what it is and the page that mentions it.
9. **More than thirty results for one study** are held for the author before they are
   accepted.

## 5. How a result is written

```
Abbreviations:
Result 1
  Vaccine:
  Who is counted:
  Number counted:
  Outcome, in the paper's words:
  Outcome, in our terms:
  By when:
  Figure:
  Other time limits reported:
Subgroups explored:
```

A result calculated from groups (rule 2) has one more line, "Calculated from:", after
"Figure".

- **Abbreviations** is one line for the study, above its results. It gives the meaning of
  every abbreviation the results use: "Abbreviations: ACIP = Advisory Committee on
  Immunization Practices; RV1 = Rotarix, a two-dose rotavirus vaccine". The results may
  then use them without explaining them again. Where there are none: "Abbreviations: none".
- **Vaccine** is the vaccine's name. Where the result is for one specific product, the
  product follows in brackets, whether or not the study compares brands: "Rotavirus
  (RV1)", "HPV (quadrivalent)", "Zoster (RZV)". Where the paper gives brands separately
  and also together: "Rotavirus (all brands)", "Rotavirus (mixed RV1 and RV5)". Where the
  result mixes products or the paper does not say, the name alone: "HPV". A vaccine reported only for comparison ends with "(comparison)": "DTaP
  (comparison)". Where one result covers several vaccines: "Combined series: polio, DTP,
  MMR" where all are required, and "Any of: DTP, Hib, IPV, PCV" where one is enough.
  Nothing else goes in this field.
- **Who is counted** begins with "everyone eligible", "first-dose recipients" or "series
  completers", then says who they are in plain words. Where the result is one group of a
  paper with no whole-population figure (rule 2), fixed words follow the opening: "largest
  group; no whole-population figure;", or, where the sizes of the groups could not be
  found, "one group; no whole-population figure; group sizes not found;". A result may refer to the group of another result of the
  same study, in the exact words "the group of Result N", followed by anything that narrows
  it: "first-dose recipients; the group of Result 1, with at least one rotavirus dose". No
  other way of pointing to another result is used.
- **Number counted** is the number of people in the group, with its page. Where no exact
  number could be found: "not found", followed by the paper's own words in quotation marks
  if it gives an approximate one: 'not found; the paper says "nearly 300,000" (p1111)'.
- **Outcome, in the paper's words** is quoted exactly, with its page. The field holds only
  the quoted words and the page. Where the paper prints only the opposite, its words for
  the opposite are quoted.
- **Outcome, in our terms** is one of: at least N doses, with the number written in ("at
  least 1 dose", "at least 2 doses"); completion; on-time completion. "Completion" is what
  the paper presents as complete, up to date or fully vaccinated, by its own definition;
  for a figure by a stated age, "By when" carries the age and the doses required: "24
  months of age (3 doses)". "At least N doses" is for a dose count the paper does not
  present as complete. A single dose of a single vaccine is always "at least 1 dose",
  whatever the paper calls it.
- **By when** is one time limit. Where no time limit could be found in the paper, it begins
  "no time limit found;" and then says how long people were observed, for example "no time
  limit found; patients followed for as long as enrolled".
- **Figure** is written as printed: numerator, denominator and percentage, with the table
  or page: "84,706 of 108,533 (78.1%) (Table 1, p4)". Where a numerator or denominator
  could not be found in the paper, the three parts are written out in order, separated by
  semicolons: "numerator not found; denominator 2,616 (p3719); 67.5% (p3719)", or
  "numerator not found; denominator not found; 32% (p1111)". Where the parts are printed
  on different pages, the pages go together at the end: "64,934 of 307,488 (21.1%)
  (numerator Table 2, p5611; denominator and percentage p5609)". The field holds one
  percentage: the paper's, where it prints an exact one; ours, marked "our calculation",
  where it prints none or only an approximate one. No quotations in this field.
- **A figure we calculated** is marked "our calculation" and names the printed numbers it
  comes from, with their pages, in this form: "119,024 of 162,614 (73.2%), our calculation:
  162,614 minus 43,590 (both Table 2, p391)".
- **A figure calculated from groups** (rule 2) is written "260,045 of 594,917 (43.7%), our
  calculation from the groups below". The line "Calculated from:" then lists each group's
  size and percentage, with the table and page: "159,916 at 50.1%; 91,352 at 39.3%; ...
  (all Table 2, pS391)". The numerator is always given: each group's size times its
  percentage, added up, and rounded once at the end.
- **Numbers** we write have commas ("31,178"), whatever the paper uses. A quotation keeps
  the paper's own form.
- **Pages** are the page number printed on the article. Where none is printed, the page of
  the PDF, written "PDF p4".
- A reader should understand the study's results without opening the paper. Every
  abbreviation is in the study's Abbreviations line, and a result does not point to
  anything outside the study's results.

### Examples

The two studies used as examples are left out of the blind verification.

**One simple result (Patrick et al. 2025)**

```
Abbreviations: RSV = respiratory syncytial virus
Result 1
  Vaccine: RSV
  Who is counted: everyone eligible; patients aged 60 and older with at
    least 1 year of membership on 23 September 2023
  Number counted: 1,003,132 (PDF p3)
  Outcome, in the paper's words: "vaccinated for RSV" (PDF p3)
  Outcome, in our terms: at least 1 dose
  By when: 9 April 2024
  Figure: 76,190 of 1,003,132 (7.6%) (Table 1, PDF p4)
  Other time limits reported: none
Subgroups explored: age, sex, race and ethnicity, Medicaid, neighbourhood
  deprivation, health conditions, prior health care use (Table 1, PDF p4)
```

**Two denominators and two brands (Calnan et al. 2016)**

```
Abbreviations: RV1 = Rotarix, a two-dose rotavirus vaccine; RV5 = RotaTeq,
  a three-dose rotavirus vaccine; ACIP = Advisory Committee on Immunization
  Practices
Result 1
  Vaccine: Rotavirus (all brands)
  Who is counted: everyone eligible; infants in four state Medicaid
    programs enrolled for eight months from birth
  Number counted: 658,219 (Table 3, p1237)
  Outcome, in the paper's words: "completed all doses" (Table 3, p1237)
  Outcome, in our terms: completion
  By when: 24 weeks of age (RV1) or 32 weeks (RV5), the package-insert
    schedule
  Figure: 178,204 of 658,219 (27.1%) (Table 3, p1237)
  Other time limits reported: none
Result 2
  Vaccine: Rotavirus (all brands)
  Who is counted: first-dose recipients; the group of Result 1, with at
    least one rotavirus dose
  Number counted: 392,957 (Table 3, p1237)
  Outcome, in the paper's words: "completed all doses" (Table 3, p1237)
  Outcome, in our terms: completion
  By when: 24 weeks of age (RV1) or 32 weeks (RV5), the package-insert
    schedule
  Figure: 178,204 of 392,957 (45.3%) (Table 3, p1237)
  Other time limits reported: none
Result 3
  Vaccine: Rotavirus (RV1)
  Who is counted: first-dose recipients; the group of Result 1, whose
    doses were RV1
  Number counted: 47,766 (Table 3, p1237)
  Outcome, in the paper's words: "completed all doses" (Table 3, p1237)
  Outcome, in our terms: completion
  By when: 24 weeks of age, the package-insert schedule
  Figure: 26,038 of 47,766 (54.5%) (Table 3, p1237)
  Other time limits reported: none
Result 4
  Vaccine: Rotavirus (RV5)
  Who is counted: first-dose recipients; the group of Result 1, whose
    doses were RV5
  Number counted: 345,191 (Table 3, p1237)
  Outcome, in the paper's words: "completed all doses" (Table 3, p1237)
  Outcome, in our terms: completion
  By when: 32 weeks of age, the package-insert schedule
  Figure: 152,166 of 345,191 (44.1%) (Table 3, p1237)
  Other time limits reported: none
Subgroups explored: year of birth, state (Table 3, p1237)
```

The study's other results follow in the same form: on-time completion ("fully compliant")
for the same four groups, and all of these again under the ACIP schedule (Table 6, p1240),
which is a different schedule and so a separate set of results.

## 6. Checks during the extraction

Done study by study, while the paper is open.

1. **Against the existing `Outcome:` line (step B).** Every change made in step B is
   listed with its reason: the new result was corrected, or the old `Outcome:` line was
   corrected. An old line that is right is left exactly as it is. A corrected line is dated
   in the study's `Verified:` line. The accuracy is documented: how many figures were
   compared, how many agreed, and, of those that differed, how many times the new results
   were right and how many times the old `Outcome:` line was right.
2. **Page check, by script.** Every figure in a result is looked for in the PDF text of the
   page the result names. A figure that is not found goes back to be read on the page.

## 7. Checks after the extraction

Done once the results of a batch are written, by agents that did not do the extraction.

1. **Consistency check.** A second agent reads every result and checks it against sections
   4 and 5. It does not open the papers. It checks that:
   - every field is present;
   - the study has its Abbreviations line;
   - "Who is counted" begins with "everyone eligible", "first-dose recipients" or "series
     completers";
   - "Outcome, in our terms" is one of the three terms;
   - the paper's words are in quotation marks and have a page, and the field holds nothing
     else;
   - the result holds one group, one outcome, one time limit and one percentage ("no time
     limit found" complies);
   - every figure has a page or table;
   - "Other time limits reported" and "Subgroups explored" hold no figures (the line
     "Calculated from:" holds the groups' sizes and percentages);
   - the wording can be understood without the paper ("the group of Result N" and an
     abbreviation given in the study's Abbreviations line comply).
   It lists every result that does not comply, and those are rewritten.
2. **Comparison with the stage 0 tables.** Another agent compares the new results with the
   stage 0 tables (`V2_P0a_Supplement_2026-10-21.docx`), for the 37 studies they cover. Every
   figure for the same quantity that differs is listed. For each, the paper is read, and
   the list records which document was right. An error in the new results is corrected.
   The accuracy is documented: how many figures were compared, how many agreed, and, of
   those that differed, how many times the new results were right and how many times the
   stage 0 tables were right.

## 8. Blind verification

- **Sample:** 30 studies drawn at random with a recorded seed, from all included studies
  except the two example studies. For a study with one result, that result is checked.
  For a study with more than one, two of its results are drawn at random. This gives
  between 30 and 60 results.
- **Checkers:** the author and a Claude Fable 5.2 agent, on the same results, neither
  seeing the other's answers.
- **What a checker sees:** the result without its figures, and the PDF. The checker enters
  the figures, or marks the result "cannot find" or "not clearly described".
- **Agreement:** a result agrees when every figure the checker enters equals the figure
  in `data_extraction.md`.
- **Threshold:** at least 80% of the results checked, for each checker. The number of
  studies with every checked result right is reported beside it.

## 9. What happens after the blind verification

- **Both at or above 80%:** the errors found are corrected and listed, the agreement of
  each checker is reported, and the stage ends. No further round.
- **Either below 80%:** stop. Find the cause, revise this procedure, and redo the
  extraction with a fresh agent.

## 10. The CSV

`study_results.csv`, generated by script from `data_extraction.md` and nothing else. One row
per result.

- **From the study's entry, as it stands in the 1D copy:** h_number, study, journal, year,
  doi_or_url, observation_requirement, denominator_type, gap_tolerance, seen_for_care,
  effect_tested, effect_discussed, exclusion, exclusion_note.
- **From the result:** vaccine_group, comparison_vaccine, who_is_counted, number_counted,
  outcome_paper_words, outcome_our_terms, by_when, numerator, denominator, percent, page,
  other_time_limits, subgroups_explored, abbreviations, largest_group_only,
  our_calculation.
- `comparison_vaccine` is "yes" where the Vaccine field ends with "(comparison)", otherwise
  "no". `largest_group_only` is "yes" where "Who is counted" says "largest group; no
  whole-population figure" or "one group; no whole-population figure", otherwise "no".
  `our_calculation` is "yes" where the Figure field says "our calculation", otherwise "no";
  the "Calculated from:" line is not carried into the CSV. `abbreviations` carries the study's Abbreviations line on each of its rows. Where a
  result says "the group of Result N", the script writes that result's description in its
  place, so that every row can be read on its own.

## 11. Order of work

1. Pilot: eleven rotavirus studies, then eleven studies of other vaccines. Extraction with
   the checks of section 6, then the checks of section 7. Report the time taken and what
   the checks found.
2. The author decides whether to go on, and whether the procedure needs changing.
3. The remaining studies, in batches of eleven, each in the same way.
4. Blind verification (section 8), once, on the finished results.
5. An independent reader checks the stage's README before the folder is sealed.

## 12. Approval and amendments

Approved by the author on 24 November 2026, and again with each amendment below after the pilots.

1. **24 November 2026, before any work began.** One agent per study (section 3).
2. **24 November 2026, after the pilot on eleven rotavirus studies.** Each brand, each
   schedule and each observation-requirement cohort is its own result, and the limit is
   twenty results (section 4, rules 3, 5, 6, 7 and 9, and the second example).
3. **24 November 2026, after the second pilot on the same eleven studies.** The example of
   a study with no result (Quinlan et al. 2015) was removed, because that study does report
   results. Added: a vaccine reported only for comparison is a result; "at least N doses";
   "no time limit found"; "numerator not found" and "denominator not found"; a result
   worked out from the printed opposite; which figure to take where the paper gives two;
   and the forms of the Vaccine field, of numbers and of calculated figures (section 4,
   rules 1 to 4, and section 5).
4. **24 November 2026, after the third pilot on the same eleven studies.** Where figures
   differ, the paper's own stated figure comes first; a numerator may be summed from
   printed counts; "the group of Result N"; the Abbreviations line; "(comparison)" in the
   Vaccine field; the form of a figure with a part not found, and of Number counted; two
   CSV columns added; a second pilot on studies of other vaccines (section 4, rules 2 and
   4; sections 5, 7, 10 and 11).
5. **24 November 2026, after the second pilot, on eleven studies of other vaccines.** A
   paper with no whole-population figure: calculated from its groups where possible,
   otherwise the largest group, marked; the data source as a reason for a separate result;
   Medicare Advantage and fee-for-service; "Combined series:" and "Any of:" in the Vaccine
   field; "completion" by the paper's own definition; "series completers"; supplements not
   in the PDF; the Figure field where parts are on different pages, and one percentage
   only; one CSV column added (section 4, rules 2, 3 and 8; sections 5, 7 and 10).
6. **24 November 2026, after the second pilot was run again.** Calculating from groups only
   from plain shares, and its form, with the line "Calculated from:"; what to take where
   group sizes could not be found; a single dose is "at least 1 dose"; an all-sources
   figure the paper states without numbers; the test for separate cohorts against several
   time limits; the product in brackets in the Vaccine field; one CSV column added
   (section 4, rules 2, 3 and 5; sections 5, 7 and 10).
