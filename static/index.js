$(document).ready(function () {
	$(".hide").click(function () {
		$(this).closest(".expandedText").hide().prev(".show").show();
		return false;
	});
	$(".show").click(function () {
		$(this).hide().next(".expandedText").show();
		return false;
	});
	$(".header").click(function () {
		$header = $(this);
		$content = $header.next();
		$content.slideToggle(100, function() {
			$header.text(function() {
				return $content.is(":visible") ? "Design primers for known alleles (collapse)" : "Design primers for known alleles";
			});
		});
	});
	$(".header2").click(function () {
		$header = $(this);
		$content = $header.next();
		$content.slideToggle(100, function() {
			$header.text(function() {
				return $content.is(":visible") ? "Design Screening Primers for Unknown CRISPR Editing (collapse)" : "Design Screening Primers for Unknown CRISPR Editing";
			});
		});
	});
});
