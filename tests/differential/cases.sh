#!/bin/bash
# Print the differential test cases, one per line:
#   suite <TAB> dataset <TAB> align_benchmark arguments (without -i/-o)
# Each case is run on Rust and C and the outputs compared (run.sh).
DATASETS="utest long"

emit() { # suite args...
  local suite=$1; shift
  for ds in $DATASETS; do printf '%s\t%s\t%s\n' "$suite" "$ds" "$*"; done
}

# Metrics and penalties; each entry is "algorithm|extra args".
METRICS=(
  "edit-wfa|" "indel-wfa|" "gap-linear-wfa|" "gap-affine-wfa|" "gap-affine2p-wfa|"
  "gap-affine-wfa|--affine-penalties=0,1,2,1"
  "gap-affine-wfa|--affine-penalties=-3,5,3,2"
)

# core: every metric x every memory mode, no heuristic (+ wfa-adaptive).
for c in "${METRICS[@]}" "gap-affine-wfa|--wfa-heuristic=wfa-adaptive --wfa-heuristic-parameters=10,50,1"; do
  for m in "" "--wfa-score-only" "--wfa-memory=med" "--wfa-memory=low" \
           "--wfa-memory=ultralow" "--wfa-memory=ultralow --wfa-score-only"; do
    emit core -a ${c%%|*} ${c#*|} $m
  done
done

# span: extension and ends-free (begin and/or end gaps), incl. match < 0.
for c in "${METRICS[@]}" "gap-affine-wfa|--affine-penalties=-5,1,2,1" \
         "gap-linear-wfa|--linear-penalties=-1,3,2" \
         "gap-affine2p-wfa|--affine2p-penalties=-1,4,6,2,24,1"; do
  for sp in extension "ends-free,0,0.1,0,0.1" "ends-free,0.1,0,0.1,0" \
            "ends-free,0.1,0.1,0.1,0.1" "ends-free,20,20,20,20"; do
    for m in "" "--wfa-score-only" "--wfa-memory=med"; do
      emit span -a ${c%%|*} ${c#*|} --wfa-span=$sp $m
    done
  done
done

# banded-static (edit/indel excluded: C never terminates, HPCBio/wfa2lib-rs#24).
for a in gap-linear-wfa gap-affine-wfa gap-affine2p-wfa; do
  for b in "-10,10" "-50,50" "-5,30"; do
    for m in "" "--wfa-score-only" "--wfa-memory=med"; do
      emit banded -a $a --wfa-heuristic=banded-static --wfa-heuristic-parameters=$b $m
    done
  done
done

# banded-adaptive
for a in gap-linear-wfa gap-affine-wfa gap-affine2p-wfa; do
  for b in "-10,10,1" "-50,50,10"; do
    for m in "" "--wfa-score-only" "--wfa-memory=med" "--wfa-memory=low"; do
      emit adaptive -a $a --wfa-heuristic=banded-adaptive --wfa-heuristic-parameters=$b $m
    done
  done
done

# drops: X-drop / Z-drop over spans and memory modes.
for h in xdrop zdrop; do
  for p in 20,1 50,1 50,10 100,5; do
    for c in "gap-affine-wfa|" "gap-affine-wfa|--affine-penalties=-3,5,3,2" \
             "gap-linear-wfa|" "gap-affine2p-wfa|"; do
      for sp in "" "--wfa-span=extension" "--wfa-span=ends-free,0.1,0.1,0.1,0.1"; do
        for m in "" "--wfa-score-only" "--wfa-memory=med"; do
          emit drops -a ${c%%|*} ${c#*|} $sp --wfa-heuristic=$h --wfa-heuristic-parameters=$p $m
        done
      done
    done
  done
done

# biwfa: BiWFA (ultralow) with every heuristic.
for h in "banded-static|-10,10" "banded-static|-50,50" "banded-adaptive|-10,10,1" \
         "banded-adaptive|-50,50,10" "xdrop|20,1" "xdrop|100,5" "zdrop|20,1" \
         "zdrop|100,5" "wfa-adaptive|10,50,1"; do
  for c in "${METRICS[@]}"; do
    a=${c%%|*}
    case "$a ${h%%|*}" in "edit-wfa banded"*|"indel-wfa banded"*|*-wfa\ xdrop|*-wfa\ zdrop)
      [ "$a" = edit-wfa ] || [ "$a" = indel-wfa ] && continue;; esac
    for m in "" "--wfa-score-only"; do
      emit biwfa -a $a ${c#*|} --wfa-memory=ultralow \
        --wfa-heuristic=${h%%|*} --wfa-heuristic-parameters=${h#*|} $m
    done
  done
done
