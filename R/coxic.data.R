### format data for coxic.c
coxic.data <- function(id, start, stop, from, to, status, z, states)
{
  p <- max(1, ncol(z))
  z <- data.frame(id, from, to, status, z)
  names(z) <- c("id", "from", "to", "status", paste("z", 1:p, sep = ""))
  ## type-specific covariates (nb: ? -> 2 presumed to hold values for 1 -> 2)
  z <- merge(merge(subset(z, from %in% states[1] & to == states[2]),
                   subset(z, from %in% states[1] & to == states[3]), by = "id"),
             subset(z, from %in% c(states[2], NA) & to == states[3]),
             by = "id", all = TRUE)
  z <- z[, substr(names(z), 1, 1) == "z"]
  z[is.na(z)] <- 0
  names(z) <- paste(paste("z", 1:p, ".", sep= ""),
                    rep(c("01", "02", "12"), each = p), sep = "")
  rownames(z) <- NULL
  uid <- unique(id)
  n <- length(uid)
  ## NA action permits missing 'start'/'stop' when 'start' = 'stop'
  start[is.na(start) & !is.na(from)] <- stop[is.na(start) & !is.na(from)]
  stop[is.na(stop)] <- start[is.na(stop)]
  ## largest and smallest observation times
  u <- as.vector(by(start, id, min, na.rm = TRUE))
  v <- as.vector(by(stop, id, max))
  ## T observed?
  absorb <- is.element(uid, id[to == states[3] & status == 1])
  ## contribution via 0 -> 1 (1), 0 -> 2 (2), both (0)?
  contrib <- rep(2, n)
  contrib[uid %in% id[from %in% states[2]]] <- 1
  contrib[uid %in% id[is.na(from)]] <- 0
  ## (possible) censoring intervals (L, R] with L = R if T01 observed exactly
  left <- stop[to == states[2]]
  left[!absorb & contrib == 2] <- v[!absorb & contrib == 2]
  right <- rep(Inf, n)
  right[contrib == 0] <- stop[is.na(from)]
  right[contrib == 1] <- start[from %in% states[2]]
  right[absorb & contrib == 0] <- v[absorb & contrib == 0] - .Machine$double.eps
  if (!any(contrib == 0) | !any(contrib == 1 & absorb))
    stop("Estimation requires some exactly-observed transition times.")
  ## maximal intersections containing 0 -> 1 support
  t01 <- maximalint(cbind(left, right)[contrib == 1, ])$int[, 2]
  names(t01) <- NULL
  t02 <- sort(unique(v[absorb & contrib == 2]))
  t12 <- sort(unique(v[absorb & contrib == 1]))
  list(supp = list(t01 = t01, t02 = t02, t12 = t12), left = left, right = right,
       u = u, v = v, contrib = contrib, absorb = absorb, z = z)
}
