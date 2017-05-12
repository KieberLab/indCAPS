$(document).ready(function () {
	$(".hide").click(function () {
		$(this).closest(".expandedText").hide().prev(".show").show();
		return false;
	});
	$(".show").click(function () {
		$(this).hide().next(".expandedText").show();
		return false;
	});
});h